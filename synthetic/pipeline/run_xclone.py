#!/usr/bin/env python3
# Run XClone (RDR + BAF) on one replicate's Standardized_Inputs.
import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite

import scanpy as sc
import xclone

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(description="Run XClone on one replicate")
    parser.add_argument("--input_dir", required=True, help="Replicate directory (Rep_*)")
    parser.add_argument("--output_dir", required=True, help="Tool output directory (Tool_Outputs/XClone)")
    parser.add_argument("--std_in_name", default="Standardized_Inputs",
                        help="Standardized input folder name under the replicate directory")
    return parser.parse_args()


def reorder_anndata_genes(adata):
    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    adata.var["chr"] = adata.var["chr"].astype(str)
    adata.var["chr"] = pd.Categorical(adata.var["chr"], categories=chrom_order, ordered=True)
    adata.var["start"] = adata.var["start"].astype(int)

    sorted_index = adata.var.sort_values(by=["chr", "start"]).index
    adata_sorted = adata[:, sorted_index].copy()
    return adata_sorted


def main():
    args = parse_args()
    input_dir = args.input_dir
    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)

    rdr_dir = os.path.join(input_dir, "synthetic", "spatial_data")
    std_xclone = os.path.join(input_dir, args.std_in_name, "XClone")
    rdr_input_dir = os.path.join(std_xclone, "RDR")
    baf_input_dir = os.path.join(std_xclone, "BAF")
    os.makedirs(rdr_input_dir, exist_ok=True)
    os.makedirs(baf_input_dir, exist_ok=True)

    anno_path = os.path.join(std_xclone, "Annotation", "annotation.csv")
    if not os.path.exists(anno_path):
        anno_path = os.path.join(std_xclone, "annotation.csv")

    # ---- RDR input preparation ----
    anno_df = pd.read_csv(anno_path)
    if "cell_id" not in anno_df.columns and "cell" in anno_df.columns:
        anno_df = anno_df.rename(columns={"cell": "cell_id"})
    anno_cells = set(anno_df["cell_id"].astype(str))

    expr_path = os.path.join(rdr_input_dir, "st-expression.parquet")
    if not os.path.exists(expr_path):
        expr_path = os.path.join(rdr_dir, "st-expression.parquet")
    expr = pd.read_parquet(expr_path)
    if "index" in expr.columns:
        expr = expr.set_index("index")  # genes x cells

    expr.columns = expr.columns.astype(str)
    expr = expr.loc[:, [c for c in expr.columns if c in anno_cells]]

    gene_anno = xclone.pp.load_anno(genome_mode="hg38_genes")
    features = gene_anno.set_index("GeneName")

    common = expr.index.intersection(features.index)
    expr = expr.loc[common]
    features = features.loc[expr.index]

    mmwrite(os.path.join(rdr_input_dir, "xcltk.rdr.mtx"), sp.csr_matrix(expr.values))
    pd.Series(expr.columns).to_csv(os.path.join(rdr_input_dir, "barcodes.tsv"), index=False, header=False)
    features.to_csv(os.path.join(rdr_input_dir, "features.tsv"), sep="\t", header=False)

    # ---- cell annotation ----
    anno_df = anno_df.rename(columns={"cell_id": "cell", "group": "cell_type"})
    barcodes = pd.read_csv(os.path.join(rdr_input_dir, "barcodes.tsv"), header=None)[0].astype(str)
    anno_df = anno_df[anno_df["cell"].astype(str).isin(set(barcodes))]
    anno_df.to_csv(os.path.join(rdr_input_dir, "cell_anno.csv"), index=False)

    # ---- Build RDR_adata ----
    RDR_adata = xclone.pp.xclonedata(
        os.path.join(rdr_input_dir, "xcltk.rdr.mtx"), data_mode="RDR",
        mtx_barcodes_file=os.path.join(rdr_input_dir, "barcodes.tsv"),
        regions_anno_file=os.path.join(rdr_input_dir, "features.tsv"),
        genome_mode="hg38_genes",
    )
    RDR_adata = xclone.pp.extra_anno(
        RDR_adata, os.path.join(rdr_input_dir, "cell_anno.csv"),
        barcodes_key="cell", cell_anno_key=["cell_type"], sep="," 
    )

    # ---- Build gene-level BAF inputs ----
    raw_baf_dir = os.path.join(std_xclone, "BAF")
    if not all(os.path.exists(os.path.join(raw_baf_dir, f)) for f in ("AD.mtx", "DP.mtx", "barcodes.tsv")):
        raw_baf_dir = os.path.join(input_dir, "baf")
    raw_AD_file = os.path.join(raw_baf_dir, "AD.mtx")
    raw_DP_file = os.path.join(raw_baf_dir, "DP.mtx")
    raw_barcodes_file = os.path.join(raw_baf_dir, "barcodes.tsv")

    full_barcodes = pd.read_csv(raw_barcodes_file, header=None)[0].astype(str)
    full_barcodes_map = pd.Series(index=full_barcodes.values, data=range(len(full_barcodes)))
    valid = full_barcodes_map.reindex(barcodes).dropna().astype(int)
    if valid.empty:
        raise ValueError("No overlapping barcodes between XClone RDR and BAF inputs.")
    valid_indices = valid.values
    barcodes = pd.Series(valid.index.astype(str))

    AD_mtx = mmread(raw_AD_file).tocsr()
    DP_mtx = mmread(raw_DP_file).tocsr()
    AD_mtx_subset = AD_mtx[valid_indices, :].T
    DP_mtx_subset = DP_mtx[valid_indices, :].T

    snp_info_path = os.path.join(raw_baf_dir, "snp_info.csv")
    if not os.path.exists(snp_info_path):
        snp_info_path = os.path.join(input_dir, "baf", "snp_info.csv")
    my_snps = pd.read_csv(snp_info_path)
    my_snps["chrom"] = my_snps["chrom"].astype(str)
    hg38_ref = xclone.pp.load_anno(genome_mode="hg38_genes")
    raw_unique_genes = my_snps["gene_id"].unique()
    ref_genes = set(hg38_ref["GeneName"].values)
    valid_genes = [g for g in raw_unique_genes if g in ref_genes]
    unique_genes = np.array(valid_genes)
    gene_to_idx = {gene: i for i, gene in enumerate(unique_genes)}

    row_indices = []
    col_indices = []
    for snp_idx, gene_name in enumerate(my_snps["gene_id"]):
        if gene_name in gene_to_idx:
            row_indices.append(snp_idx)
            col_indices.append(gene_to_idx[gene_name])

    agg_operator = sp.coo_matrix(
        (np.ones(len(row_indices)), (row_indices, col_indices)),
        shape=(len(my_snps), len(unique_genes))
    ).tocsr()

    AD_gene = AD_mtx_subset.T @ agg_operator
    DP_gene = DP_mtx_subset.T @ agg_operator
    final_gene_features = hg38_ref.set_index("GeneName").reindex(unique_genes).reset_index()

    mmwrite(os.path.join(baf_input_dir, "AD.mtx"), AD_gene.T)
    mmwrite(os.path.join(baf_input_dir, "DP.mtx"), DP_gene.T)
    barcodes.to_csv(os.path.join(baf_input_dir, "barcodes.tsv"), index=False, header=False)
    final_gene_features.to_csv(os.path.join(baf_input_dir, "features.tsv"), sep="\t", index=False, header=False)

    # ---- Build BAF_adata ----
    BAF_adata = xclone.pp.xclonedata(
        [os.path.join(baf_input_dir, "AD.mtx"), os.path.join(baf_input_dir, "DP.mtx")],
        data_mode="BAF",
        mtx_barcodes_file=os.path.join(baf_input_dir, "barcodes.tsv"),
        regions_anno_file=os.path.join(baf_input_dir, "features.tsv"),
        genome_mode="hg38_genes",
    )
    BAF_adata = xclone.pp.extra_anno(
        BAF_adata, os.path.join(rdr_input_dir, "cell_anno.csv"),
        barcodes_key="cell", cell_anno_key=["cell_type"], sep="," 
    )

    # ---- Align RDR/BAF gene metadata and ordering ----
    column_names = ["GeneID", "chr", "start", "stop", "arm", "chr_arm", "band"]
    RDR_adata.var.columns = column_names
    RDR_adata.var.index.name = "GeneName"
    RDR_adata.var["GeneName"] = RDR_adata.var.index

    BAF_adata.var.columns = column_names
    BAF_adata.var.index.name = "GeneName"
    BAF_adata.var["GeneName"] = BAF_adata.var.index

    baf_genes = set(BAF_adata.var_names)
    rdr_genes = set(RDR_adata.var_names)
    common_genes = list(baf_genes.intersection(rdr_genes))

    temp_var = RDR_adata.var.loc[common_genes]
    sorted_common_genes = temp_var.sort_values(by=["chr", "start"]).index

    BAF_adata = BAF_adata[:, sorted_common_genes].copy()
    RDR_adata = RDR_adata[:, sorted_common_genes].copy()

    RDR_adata = reorder_anndata_genes(RDR_adata)
    BAF_adata = reorder_anndata_genes(BAF_adata)

    RDR_adata.var.reset_index(drop=True, inplace=True)
    BAF_adata.var.reset_index(drop=True, inplace=True)

    # ---- Run XClone ----
    dataset_name = os.path.basename(input_dir)
    xclone.pp.efficiency_preview()

    xc_rdr = xclone.XCloneConfig(dataset_name=dataset_name, module="RDR")
    xc_rdr.outdir = outdir
    xc_rdr.cell_anno_key = "cell_type"
    xc_rdr.ref_celltype = "Normal"
    xc_rdr.marker_group_anno_key = "cell_type"
    xc_rdr.top_n_marker = 25
    xc_rdr.xclone_plot = False
    xc_rdr.display()

    RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file=xc_rdr)

    xc_baf = xclone.XCloneConfig(dataset_name=dataset_name, module="BAF")
    xc_baf.outdir = outdir
    xc_baf.cell_anno_key = "cell_type"
    xc_baf.ref_celltype = "Normal"
    xc_baf.update_info_from_rdr = False
    xc_baf.BAF_denoise = True
    xc_baf.xclone_plot = False
    xc_baf.display()

    BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, config_file=xc_baf)

    xc_comb = xclone.XCloneConfig(dataset_name=dataset_name, module="Combine")
    xc_comb.outdir = outdir
    xc_comb.cell_anno_key = "cell_type"
    xc_comb.ref_celltype = "Normal"
    xc_comb.BAF_denoise = True
    xc_comb.xclone_plot = True
    xc_comb.display()

    combine_Xdata = xclone.model.run_combine(
        RDR_Xdata, BAF_merge_Xdata, verbose=True, run_verbose=True, config_file=xc_comb
    )

    # ---- Export CNV state matrix ----
    combine_path = os.path.join(outdir, "data", "combined_final.h5ad")
    combine_Xdata = sc.read(combine_path)

    use_layer = "prob1_merge"  # or "corrected_prob"
    X = combine_Xdata.layers[use_layer]

    if X.ndim == 4:
        X_rdr = X.sum(axis=-1)
        cna_state = X_rdr.argmax(axis=2)
    elif X.ndim == 3:
        cna_state = X.argmax(axis=2)
    else:
        raise ValueError(f"Unexpected ndim={X.ndim} for layer {use_layer}")

    state_map = {0: "loss", 1: "loh", 2: "neutral", 3: "gain"}
    state_labels = np.vectorize(state_map.get)(cna_state)

    genes = combine_Xdata.var["GeneName"] if "GeneName" in combine_Xdata.var.columns else combine_Xdata.var_names
    df = pd.DataFrame(state_labels, index=combine_Xdata.obs_names, columns=genes)

    out_tsv = os.path.join(outdir, f"{use_layer}_cnv_state_spot_by_gene.tsv")
    df.to_csv(out_tsv, sep="\t")

    print(f"[Output] CNV state matrix saved: {out_tsv}")
    print("[Done] XClone pipeline finished.")


if __name__ == "__main__":
    main()
