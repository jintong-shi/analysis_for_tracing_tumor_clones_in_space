#!/usr/bin/env python3
"""infercnvpy CLI wrapper for simulated (or real) ST data.

Takes either an .h5ad or a gene-by-cell expression Parquet + annotation Parquet,
defines the reference from a spot-ID list or an annotation group label, merges a
gene order file, normalizes (CPM + log1p), and writes the CNV matrix to Parquet.

  python run_infercnvpy.py --data_dir <sim_run> --gene_order_file <gene_pos.txt> \
      --ref_group Benign --out_dir <sim_run>/infercnvpy_out
"""

import os
import sys
import csv
import json
import argparse
from datetime import datetime
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

# Ensure writable temp/cache dirs for restricted environments
_tmp = os.environ.get("TMPDIR", "/tmp")
os.environ.setdefault("TMPDIR", _tmp)
os.environ.setdefault("TEMP", _tmp)
os.environ.setdefault("TMP", _tmp)
os.environ.setdefault("JOBLIB_TEMP_FOLDER", _tmp)

# ---------------------------
# Utilities
# ---------------------------

def log(msg: str):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] {msg}")


def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)


def load_gene_order_table(fp: str, gene_col: str = "gene", chr_col: str = "chromosome",
                          start_col: str = "start", end_col: str = "end") -> pd.DataFrame:
    """Load gene order file and return a DataFrame indexed by gene symbol/id."""
    log(f"Reading gene order from: {fp}")
    go = pd.read_table(
        fp,
        header=None if gene_col == "gene" and chr_col == "chromosome" else 0,
        quoting=csv.QUOTE_MINIMAL,
        on_bad_lines="warn",
        dtype={start_col: 'Int64', end_col: 'Int64'}
    )

    # If no header provided: set default names
    if go.shape[1] >= 4 and (go.columns.dtype == 'int64' or list(go.columns) == list(range(go.shape[1]))):
        go = go.iloc[:, :4]
        go.columns = [gene_col, chr_col, start_col, end_col]

    # Keep only required columns and drop NA genes
    cols_needed = [gene_col, chr_col, start_col, end_col]
    missing = [c for c in cols_needed if c not in go.columns]
    if missing:
        raise ValueError(f"Gene order is missing columns: {missing}")
    go = go[cols_needed].dropna(subset=[gene_col])

    # Harmonize chromosome labels
    go[chr_col] = go[chr_col].astype(str).str.replace('^chr', '', regex=True)
    go[chr_col] = np.where(go[chr_col].str.upper().isin(['M', 'MT']), 'chrM', 'chr' + go[chr_col])

    go = go.drop_duplicates(subset=[gene_col]).set_index(gene_col)
    return go


def make_anndata_from_parquet(data_dir: str, anno_fp: str | None = None) -> sc.AnnData:
    """Build AnnData from gene-by-cell parquet (rows=genes, columns=cells),
    matching the layout used in your R script.
    - If anno_fp is None, attempt to read `<data_dir>/spatial_data/st-annotation.parquet`.
    """
    expr_fp = os.path.join(data_dir, "spatial_data", "st-expression.parquet")
    if not os.path.exists(expr_fp):
        raise FileNotFoundError(f"Expression parquet not found: {expr_fp}")

    log(f"Reading expression parquet: {expr_fp}")
    mat = pd.read_parquet(expr_fp)
    if 'index' in mat.columns:
        mat = mat.set_index('index')

    if anno_fp is None:
        anno_fp = os.path.join(data_dir, "spatial_data", "st-annotation.parquet")
    if not os.path.exists(anno_fp):
        raise FileNotFoundError(f"Annotation parquet not found: {anno_fp}")

    log(f"Reading annotation parquet: {anno_fp}")
    anno = pd.read_parquet(anno_fp)

    # Standardize index name for cell IDs
    if 'cell_id' in anno.columns:
        anno = anno.set_index('cell_id')
    elif 'index' in anno.columns:
        anno = anno.set_index('index')

    # Reorder matrix columns to match annotation cells
    common_cells = [c for c in anno.index if c in mat.columns]
    if len(common_cells) == 0:
        raise ValueError("No overlapping cell/spot IDs between expression and annotation.")

    if len(common_cells) < len(anno.index):
        warnings.warn(f"Annotation has {len(anno.index)} cells, but only {len(common_cells)} found in matrix. Subsetting to overlap.")

    mat = mat.loc[:, common_cells]

    # AnnData expects cells x genes
    adata = sc.AnnData(X=mat.T.values)
    adata.var_names = mat.index.astype(str)
    adata.obs_names = list(mat.columns)

    # Bring all annotation columns into obs
    anno = anno.loc[common_cells]
    for col in anno.columns:
        adata.obs[col] = anno[col].astype(str)

    return adata


def attach_gene_order(adata: sc.AnnData, gene_order: pd.DataFrame,
                      chr_col: str = "chromosome",
                      start_col: str = "start", end_col: str = "end") -> sc.AnnData:
    """Subset to common genes and attach genomic coordinates into adata.var."""
    common = gene_order.index.intersection(adata.var_names)
    if len(common) < 2000:
        warnings.warn(f"Only {len(common)} common genes between adata and gene_order; results may be unstable.")
    adata = adata[:, list(common)].copy()

    go_sub = gene_order.loc[common, [chr_col, start_col, end_col]].copy()
    go_sub.columns = ["chromosome", "start", "end"]  # enforce names infercnvpy expects

    var = pd.DataFrame(index=adata.var_names)
    adata.var = var.merge(go_sub, left_index=True, right_index=True, how='left')

    # Drop genes missing coordinates
    miss = adata.var['chromosome'].isna().sum()
    if miss > 0:
        warnings.warn(f"{miss} genes missing coordinates were dropped.")
        keep = adata.var['chromosome'].notna()
        adata = adata[:, keep].copy()

    return adata


def add_reference_labels(adata: sc.AnnData, ref_spots: list[str] | None,
                         anno_ref_col: str | None, ref_group: str | None,
                         label_key: str = "ref_label"):
    """Create a binary reference label column in adata.obs."""
    if ref_spots is not None:
        ref_spots = set(map(str, ref_spots))
        adata.obs[label_key] = ["Ref" if s in ref_spots else "other" for s in adata.obs_names]
        return label_key, ["Ref"]

    if anno_ref_col and ref_group:
        if anno_ref_col not in adata.obs.columns:
            raise ValueError(f"Annotation column '{anno_ref_col}' not in adata.obs.")
        adata.obs[label_key] = np.where(adata.obs[anno_ref_col].astype(str) == str(ref_group), "Ref", "other")
        return label_key, ["Ref"]

    # No reference: create a dummy label so infercnvpy can still run in no-ref mode by passing empty list
    label_key_fallback = label_key
    adata.obs[label_key_fallback] = "all"
    return label_key_fallback, []


# ---------------------------
# Main CLI
# ---------------------------

def main():
    p = argparse.ArgumentParser(description="Run infercnvpy reproducibly on ST (sim or real) data.")

    io = p.add_argument_group("Inputs")
    io.add_argument("--input_h5ad", type=str, default=None, help="Existing AnnData (.h5ad). Alternative to --data_dir mode.")
    io.add_argument("--data_dir", type=str, default=None, help="Root of sim data. Expects spatial_data/st-expression.parquet and st-annotation.parquet unless --ref_annotation_file provided.")
    io.add_argument("--ref_annotation_file", type=str, default=None, help="Optional Parquet with columns 'cell_id' (or 'index') and a reference group column.")
    io.add_argument("--ref_spots_file", type=str, default=None, help="Optional plain-text file with one reference spot ID per line.")

    go = p.add_argument_group("Gene order")
    go.add_argument("--gene_order_file", type=str, required=True, help="Tab file with columns: gene, chromosome, start, end (or provide names via --gene_col etc.)")
    go.add_argument("--gene_col", type=str, default="gene")
    go.add_argument("--chr_col", type=str, default="chromosome")
    go.add_argument("--start_col", type=str, default="start")
    go.add_argument("--end_col", type=str, default="end")

    ref = p.add_argument_group("Reference definition")
    ref.add_argument("--anno_ref_col", type=str, default=None, help="Column in annotation that holds group labels (e.g. 'cell_type').")
    ref.add_argument("--ref_group", type=str, default=None, help="Name of the benign/normal reference group in the annotation column.")

    out = p.add_argument_group("Outputs")
    out.add_argument("--out_dir", type=str, required=True)
    out.add_argument("--out_prefix", type=str, default="infercnvpy")
    #out.add_argument("--save_heatmap", action="store_true", help="If set, generate chromosome heatmap PNG.")

    prep = p.add_argument_group("Preprocessing")
    prep.add_argument("--target_sum", type=float, default=1e6)
    prep.add_argument("--no_log1p", action="store_true")

    run = p.add_argument_group("infercnvpy run")
    run.add_argument("--exclude_chromosomes", type=str, nargs='*', default=["chrM"]) 
    run.add_argument("--window_size", type=int, default=100)
    run.add_argument("--n_jobs", type=int, default=8)

    args = p.parse_args()

    if (args.input_h5ad is None) == (args.data_dir is None):
        p.error("Provide exactly one of --input_h5ad or --data_dir")

    ensure_outdir(args.out_dir)

    # 1) Build AnnData
    if args.input_h5ad:
        log(f"Reading h5ad: {args.input_h5ad}")
        adata = sc.read_h5ad(args.input_h5ad)
    else:
        adata = make_anndata_from_parquet(args.data_dir, args.ref_annotation_file)

    log(f"Initial AnnData: cells={adata.n_obs}, genes={adata.n_vars}")

    # 2) Reference labels
    ref_spots = None
    if args.ref_spots_file and os.path.exists(args.ref_spots_file):
        with open(args.ref_spots_file) as f:
            ref_spots = [ln.strip() for ln in f if ln.strip()]
        log(f"Loaded {len(ref_spots)} reference spot IDs from {args.ref_spots_file}")

    ref_key, ref_cat = add_reference_labels(
        adata,
        ref_spots=ref_spots,
        anno_ref_col=args.anno_ref_col,
        ref_group=args.ref_group,
        label_key="ref_label",
    )
    if ref_cat:
        counts = adata.obs[ref_key].value_counts().to_dict()
        log(f"Reference summary: {counts}")
    else:
        log("Running in no-reference mode (ref_cat=[]).")

    # 3) Gene order & attach
    gene_order = load_gene_order_table(
        args.gene_order_file,
        gene_col=args.gene_col,
        chr_col=args.chr_col,
        start_col=args.start_col,
        end_col=args.end_col,
    )
    adata = attach_gene_order(adata, gene_order)
    log(f"After attaching gene order: cells={adata.n_obs}, genes={adata.n_vars}")

    # 4) Normalize & log1p
    sc.pp.normalize_total(adata, target_sum=float(args.target_sum))
    if not args.no_log1p:
        sc.pp.log1p(adata)
    log("Normalization complete.")

    # Save input subset for reproducibility
    input_sub_h5ad = os.path.join(args.out_dir, f"{args.out_prefix}_input.h5ad")
    adata.write(input_sub_h5ad)
    log(f"Saved subsetted input → {input_sub_h5ad}")

    # 5) Run infercnvpy
    log("Running infercnvpy...")
    ref_key_pass = None if not ref_cat else ref_key
    ref_cat_pass = None if not ref_cat else ref_cat

    cnv.tl.infercnv(
        adata,
        reference_key=ref_key_pass,
        reference_cat=ref_cat_pass,
        #reference_key=ref_key,
        #reference_cat=ref_cat,  # [] means no-ref mode
        exclude_chromosomes=args.exclude_chromosomes,
        calculate_gene_values=True,
        window_size=int(args.window_size),
        n_jobs=int(args.n_jobs),
    )
    log("infercnvpy complete.")

    # 6) Save outputs
    # out_h5ad = os.path.join(args.out_dir, f"{args.out_prefix}.h5ad")
    # adata.write(out_h5ad)
    # log(f"Saved inferCNV result → {out_h5ad}")

    # # Gene-level CNV matrix
    # cnv_mat = None
    # if "gene_values_cnv" in adata.layers:
    #     cnv_mat = adata.layers["gene_values_cnv"]  # genes x cells after transpose below
    #     df = pd.DataFrame(cnv_mat.T, index=adata.var_names, columns=adata.obs_names)
    # elif "X_cnv" in adata.obsm:
    #     # Some versions store in obsm, shape is cells x genes
    #     mat = adata.obsm["X_cnv"]
    #     df = pd.DataFrame(mat.T, index=adata.var_names, columns=adata.obs_names)
    # else:
    #     raise RuntimeError("CNV matrix not found in adata.layers['gene_values_cnv'] or adata.obsm['X_cnv'].")

    # parquet_fp = os.path.join(args.out_dir, f"{args.out_prefix}_gene_cnv.parquet")
    # df.reset_index(names="gene").to_parquet(parquet_fp, compression="snappy")
    # log(f"Saved gene-level CNV matrix → {parquet_fp}")

    # 6) Save outputs
    out_h5ad = os.path.join(args.out_dir, f"{args.out_prefix}.h5ad")
    adata.write(out_h5ad)
    log(f"Saved inferCNV result → {out_h5ad}")

    # Gene-level CNV matrix (robust)
    cnv_df = None
    try:
        if "gene_values_cnv" in adata.layers and adata.layers["gene_values_cnv"] is not None:
            mat = adata.layers["gene_values_cnv"]  # cells x genes or genes x cells, depending on version
            mat = np.asarray(mat)
            # infercnvpy usually gives cells x genes; we want genes x cells
            if mat.shape[0] == adata.n_obs and mat.shape[1] == adata.n_vars:
                cnv_df = pd.DataFrame(mat.T, index=adata.var_names, columns=adata.obs_names)
            elif mat.shape[0] == adata.n_vars and mat.shape[1] == adata.n_obs:
                cnv_df = pd.DataFrame(mat, index=adata.var_names, columns=adata.obs_names)
            else:
                raise ValueError(f"Unexpected CNV shape from layers['gene_values_cnv']: {mat.shape}")
        elif "X_cnv" in adata.obsm and adata.obsm["X_cnv"] is not None:
            mat = np.asarray(adata.obsm["X_cnv"])  # cells x genes
            if mat.shape[0] != adata.n_obs or mat.shape[1] != adata.n_vars:
                raise ValueError(f"Unexpected CNV shape from obsm['X_cnv']: {mat.shape}")
            cnv_df = pd.DataFrame(mat.T, index=adata.var_names, columns=adata.obs_names)
        else:
            warnings.warn("CNV matrix not found in adata.layers['gene_values_cnv'] nor adata.obsm['X_cnv'].")
    except Exception as e:
        warnings.warn(f"Failed to extract CNV matrix: {e}")

    if cnv_df is not None:
        parquet_fp = os.path.join(args.out_dir, f"{args.out_prefix}_gene_cnv.parquet")
        cnv_df.reset_index(names="gene").to_parquet(parquet_fp, compression="snappy")
        log(f"Saved gene-level CNV matrix → {parquet_fp}")

    # Optional heatmap
    #try:
    groupby_arg = (ref_key if ref_cat else None)
    fig = cnv.pl.chromosome_heatmap(adata, groupby=groupby_arg, dendrogram=True)
    heatmap_fp = os.path.join(args.out_dir, f"{args.out_prefix}_heatmap.png")
    plt.savefig(heatmap_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved heatmap → {heatmap_fp}")
    
    # Write a small run metadata JSON
    meta = {
        "input_mode": "h5ad" if args.input_h5ad else "parquet",
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "exclude_chromosomes": args.exclude_chromosomes,
        "window_size": int(args.window_size),
        "n_jobs": int(args.n_jobs),
        "target_sum": float(args.target_sum),
        "log1p": not args.no_log1p,
        "ref_counts": (adata.obs[ref_key].value_counts().to_dict() if (ref_key and ref_key in adata.obs) else {}),
    }
    with open(os.path.join(args.out_dir, f"{args.out_prefix}_runmeta.json"), 'w') as f:
        json.dump(meta, f, indent=2)
    log("Run metadata written.")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        log(f"ERROR: {e}")
        sys.exit(1)
