#!/usr/bin/env python3
# CopyVAE input prep: counts.mtx + genes.tsv + barcodes.tsv -> st_expression.h5ad
import scanpy as sc
import pandas as pd
import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Directory containing counts.mtx, genes.tsv, barcodes.tsv")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_h5ad = os.path.join(input_dir, "st_expression.h5ad")
    
    files = ["counts.mtx", "genes.tsv", "barcodes.tsv"]
    for f in files:
        if not os.path.exists(os.path.join(input_dir, f)):
            print(f"Error: Missing {f} in {input_dir}")
            sys.exit(1)

    print(f"[Prep] Converting MTX to H5AD in {input_dir}...")

    # writeMM exports Genes x Cells; scanpy needs Cells x Genes
    adata = sc.read_mtx(os.path.join(input_dir, "counts.mtx")).T 
    
    genes = pd.read_csv(os.path.join(input_dir, "genes.tsv"), header=None)[0].values
    barcodes = pd.read_csv(os.path.join(input_dir, "barcodes.tsv"), header=None)[0].values

    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()

    # CopyVAE looks up var['gene_name'] and errors out if only the index is set
    adata.var['gene_name'] = adata.var_names

    adata.write(output_h5ad)
    print(f"[Prep] Saved: {output_h5ad}")

if __name__ == "__main__":
    main()