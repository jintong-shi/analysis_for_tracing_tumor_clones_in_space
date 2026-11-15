#!/usr/bin/env python3
# CopyVAE post-process: copy.npy -> log2 spots x genes CSV.
import argparse
import os
import sys
import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description="Process CopyVAE outputs into log2 spots x genes CSV.")
    ap.add_argument("--tool_dir", help="Tool_Outputs/CopyVAE directory")
    ap.add_argument("--batch_path", help="Rep_* directory (will use Tool_Outputs/CopyVAE under it)")
    ap.add_argument("--eps", type=float, default=1e-12, help="Small epsilon for log2")
    ap.add_argument("--out_file", default=None, help="Optional override output CSV")
    return ap.parse_args()


def main():
    args = parse_args()
    tool_dir = args.tool_dir
    if tool_dir is None:
        if args.batch_path is None:
            print("ERROR: Provide --tool_dir or --batch_path", file=sys.stderr)
            sys.exit(2)
        tool_dir = os.path.join(args.batch_path, "Tool_Outputs", "CopyVAE")

    if not os.path.isdir(tool_dir):
        print(f"[Skip] CopyVAE dir not found: {tool_dir}")
        sys.exit(0)

    f_copy = os.path.join(tool_dir, "copy.npy")
    f_map = os.path.join(tool_dir, "bin_to_gene_map.csv")
    f_barc = os.path.join(tool_dir, "barcodes_order.txt")
    out_file = args.out_file or os.path.join(tool_dir, "processed_log2_spots_x_genes.csv")

    if not (os.path.exists(f_copy) and os.path.exists(f_map)):
        print(f"[Skip] Missing copy.npy or bin_to_gene_map.csv in {tool_dir}")
        sys.exit(0)

    copy_mat = np.load(f_copy)  # cells x bins
    if copy_mat.ndim != 2:
        raise ValueError(f"copy.npy must be 2D, got shape {copy_mat.shape}")

    barcodes = None
    if os.path.exists(f_barc):
        with open(f_barc, "r", encoding="utf-8") as f:
            barcodes = [line.strip() for line in f if line.strip()]
        if len(barcodes) != copy_mat.shape[0]:
            print(f"[Warn] barcodes_order.txt length ({len(barcodes)}) != nrows ({copy_mat.shape[0]}). Ignoring barcodes.")
            barcodes = None

    # bins x cells, baseline 2 -> ratio
    bins_x_cells = (copy_mat / 2.0).T

    mapping = pd.read_csv(f_map)
    if not {"gene", "bin_index"}.issubset(mapping.columns):
        raise ValueError("bin_to_gene_map.csv must contain columns: gene, bin_index")

    ridx = mapping["bin_index"].astype(int).to_numpy()
    nbins = bins_x_cells.shape[0]
    ok = (ridx >= 0) & (ridx < nbins)
    if not ok.all():
        mapping = mapping.loc[ok].reset_index(drop=True)
        ridx = mapping["bin_index"].astype(int).to_numpy()
        print(f"[Warn] Dropped {np.sum(~ok)} out-of-range bin_index rows")

    gene_x_cells = bins_x_cells[ridx, :]
    gene_x_cells = pd.DataFrame(gene_x_cells, index=mapping["gene"].astype(str))

    if gene_x_cells.index.duplicated().any():
        gene_x_cells = gene_x_cells.groupby(level=0).mean()

    # log2 ratio
    log2_mat = np.log2(gene_x_cells + args.eps)

    # Spots x Genes
    final_mat = log2_mat.T
    if barcodes is not None and len(barcodes) == final_mat.shape[0]:
        final_mat.index = barcodes

    final_mat.to_csv(out_file, index=True)
    print(f"[CopyVAE] Saved: {out_file} (spots x genes = {final_mat.shape})")


if __name__ == "__main__":
    main()
