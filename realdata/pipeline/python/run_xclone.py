import sys
import os
import argparse
import pandas as pd
import numpy as np
import scipy
import scanpy as sc
import xclone
import anndata as ad

def parse_args():
    parser = argparse.ArgumentParser(description="Run XClone Analysis Pipeline")
    # Input Directories from xcltk output
    parser.add_argument("--rdr_dir", required=True, help="Path to RDR output directory containing matrix.mtx")
    parser.add_argument("--baf_dir", required=True, help="Path to BAF output directory (e.g. .../BAF/3_baf_fc/)")
    
    # Barcode Lists
    parser.add_argument("--keep_barcodes", required=True, help="Path to keep_barcodes.txt (Tumor)")
    parser.add_argument("--ref_barcodes", required=True, help="Path to ref_barcodes.txt (Normal)")
    
    # Output and Config
    parser.add_argument("--outdir", required=True, help="Final output directory")
    parser.add_argument("--dataset_name", default="XClone_Analysis")
    parser.add_argument("--genome_mode", default="hg38_genes", help="e.g., hg38_genes or hg19_genes")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    
    # ==========================================
    # 1. Load RDR Data
    # ==========================================
    print("[XClone] Loading RDR Data...")
    rdr_file = os.path.join(args.rdr_dir, "matrix.mtx")
    rdr_barcodes = os.path.join(args.rdr_dir, "barcodes.tsv")
    
    # Note: xcltk output usually needs features.tsv. If xclonedata requires it in the dir:
    # We assume features.tsv is present in rdr_dir based on standard xcltk output
    RDR_adata = xclone.pp.xclonedata(rdr_file, 'RDR', rdr_barcodes, genome_mode=args.genome_mode)
    
    # ==========================================
    # 2. Preprocessing & Connectivities
    # ==========================================
    print(f"[XClone] Preprocessing RDR. Cells before filter: {RDR_adata.n_obs}")
    
    # Standard Scanpy filtering
    sc.pp.filter_genes(RDR_adata, min_cells=1)
    sc.pp.filter_cells(RDR_adata, min_counts=100)
    print(f"[XClone] Cells after filter: {RDR_adata.n_obs}")

    # Calculate Neighbors (Required for XClone Smoothing)
    # We use a copy for calculation to avoid altering raw counts too much before XClone
    sc.pp.normalize_total(RDR_adata)
    sc.pp.log1p(RDR_adata)
    sc.pp.pca(RDR_adata)
    sc.pp.neighbors(RDR_adata)    

    # ==========================================
    # 3. Load BAF Data
    # ==========================================
    print("[XClone] Loading BAF Data...")
    ad_file = os.path.join(args.baf_dir, "xcltk.AD.mtx")
    dp_file = os.path.join(args.baf_dir, "xcltk.DP.mtx")
    baf_barcodes = os.path.join(args.baf_dir, "xcltk.samples.tsv")
    
    BAF_adata = xclone.pp.xclonedata([ad_file, dp_file], 'BAF', baf_barcodes, genome_mode=args.genome_mode)

    # ==========================================
    # 4. Annotate Cells (Normal vs Tumor)
    # ==========================================
    print("[XClone] Annotating cells...")
    normal_df = pd.read_csv(args.ref_barcodes, header=None, names=["cell"])
    tumor_df = pd.read_csv(args.keep_barcodes, header=None, names=["cell"])
    
    normal_df["cluster.pred"] = "N"
    tumor_df["cluster.pred"] = "T"
    
    anno_df = pd.concat([normal_df, tumor_df], ignore_index=True).set_index("cell")
    
    # Merge annotation
    RDR_adata.obs = RDR_adata.obs.join(anno_df, how="left")
    BAF_adata.obs = BAF_adata.obs.join(anno_df, how="left")
    
    # Filter to only keep cells present in the annotation lists (Intersection)
    # This aligns RDR, BAF, and the specific barcodes you want to analyze
    valid_cells = anno_df.index.intersection(RDR_adata.obs_names).intersection(BAF_adata.obs_names)
    RDR_adata = RDR_adata[valid_cells, :].copy()
    BAF_adata = BAF_adata[valid_cells, :].copy()
    print(f"[XClone] Final aligned cell count: {len(valid_cells)}")

    # ==========================================
    # 5. Run Modules
    # ==========================================
    
    # --- RDR Module ---
    print("[XClone] Running RDR Module...")
    xconfig_rdr = xclone.XCloneConfig(dataset_name=args.dataset_name, module="RDR")
    xconfig_rdr.outdir = args.outdir
    xconfig_rdr.cell_anno_key = "cluster.pred"
    xconfig_rdr.ref_celltype = "N"
    xconfig_rdr.xclone_plot = False # Disable plotting display
    # xconfig_rdr.display() # Disabled
    
    RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file=xconfig_rdr)
    
    # --- BAF Module ---
    print("[XClone] Running BAF Module...")
    xconfig_baf = xclone.XCloneConfig(dataset_name=args.dataset_name, module="BAF")
    xconfig_baf.outdir = args.outdir
    xconfig_baf.cell_anno_key = "cluster.pred"
    xconfig_baf.ref_celltype = "N"
    xconfig_baf.xclone_plot = False
    
    # Important: Pass RDR_adata to use its connectivities for smoothing
    BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, config_file=xconfig_baf)
    
    # --- Combine Module ---
    print("[XClone] Running Combine Module...")
    xconfig_comb = xclone.XCloneConfig(dataset_name=args.dataset_name, module="Combine")
    xconfig_comb.outdir = args.outdir
    xconfig_comb.cell_anno_key = "cluster.pred"
    xconfig_comb.ref_celltype = "N"
    xconfig_comb.xclone_plot = False
    xconfig_comb.BAF_denoise = True
    
    # Intersect chr sets: BAF retains chrY from cellsnp pileup but RDR often drops it
    # because no chrY gene survives the expression filter, which trips XClone's strict
    # chr-equality check in check_RDR_BAF_chrmapping. Dropping the orphan chr is safe
    # (gene-level chrY CNV signal is unreliable on Visium anyway).
    common_chrs = set(RDR_Xdata.var["chr"]) & set(BAF_merge_Xdata.var["chr"])
    rdr_only = set(RDR_Xdata.var["chr"]) - common_chrs
    baf_only = set(BAF_merge_Xdata.var["chr"]) - common_chrs
    if rdr_only or baf_only:
        print(f"[XClone] dropping non-shared chrs (RDR-only={sorted(rdr_only)}, BAF-only={sorted(baf_only)})")
        RDR_Xdata       = RDR_Xdata[:,       RDR_Xdata.var["chr"].isin(common_chrs)].copy()
        BAF_merge_Xdata = BAF_merge_Xdata[:, BAF_merge_Xdata.var["chr"].isin(common_chrs)].copy()

    xclone.model.run_combine(RDR_Xdata, BAF_merge_Xdata, verbose=True, config_file=xconfig_comb)
    
    # Reload the combined data
    combined_h5ad_path = os.path.join(args.outdir, "data", "combined_final.h5ad")
    merge_data = sc.read_h5ad(combined_h5ad_path)

    # ==========================================
    # 6. Save Results (Labeled Matrix & Probs)
    # ==========================================
    print("[XClone] Saving Output Matrices...")
    
    # A. Save Labeled CNV Matrix
    # 0, 1, 2, 3 -> Labeled
    cnv_states_matrix = np.argmax(merge_data.layers['prob1_merge'], axis=2)
    
    df_data = pd.DataFrame(
        cnv_states_matrix.T, 
        index=merge_data.var.index, 
        columns=merge_data.obs_names
    )
    
    # Mapping based on your provided logic
    state_map = {0: "Loss", 1: "LOH", 2: "Neutral", 3: "Gain"}
    df_data = df_data.replace(state_map)
    
    df_meta = merge_data.var[['GeneName', 'chr', 'start', 'stop']]
    final_df = pd.concat([df_meta, df_data], axis=1).set_index('GeneName')
    
    outfile_labeled = os.path.join(args.outdir, "CNV_states_matrix_with_coords_labeled.csv")
    final_df.to_csv(outfile_labeled)
    print(f"Saved Labeled Matrix to: {outfile_labeled}")

    # B. Save Probability Matrices for each state
    # Using extract_xclone_matrix as in your script
    use_layer = "prob1_merge"
    region_lst = ["GeneName", "GeneID"]
    
    states_dict = {
        0: "prob_combine_copyloss",
        1: "prob_combine_loh",
        2: "prob_combine_copyneutral",
        3: "prob_combine_copygain"
    }
    
    for state_idx, folder_name in states_dict.items():
        out_sub = os.path.join(args.outdir, folder_name + "/")
        print(f"Saving probability for state {state_idx} ({folder_name})...")
        xclone.al.extract_xclone_matrix(
            merge_data, 
            use_layer, 
            states=state_idx, 
            region_lst=region_lst, 
            out_dir=out_sub
        )

    # Create DONE flag
    with open(os.path.join(args.outdir, "_DONE.txt"), "w") as f:
        f.write("XClone Pipeline Completed Successfully.")

if __name__ == "__main__":
    main()