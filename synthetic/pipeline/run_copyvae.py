#!/usr/bin/env python3
# Run CopyVAE on st_expression.h5ad.
import os
import sys
import subprocess
import argparse

CELL_CYCLE_FILE = "<TOOLS_ROOT>/copyVAE/data/Macosko_cell_cycle_genes.txt"
GENE_MAP_CSV = "<REF_ROOT>/hg38_gene_map.csv"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Directory containing st_expression.h5ad")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    h5ad_file = os.path.join(args.input_dir, "st_expression.h5ad")
    
    if not os.path.exists(h5ad_file):
        print(f"Error: Input file missing: {h5ad_file}")
        sys.exit(1)
        
    if not os.path.exists(CELL_CYCLE_FILE):
        print(f"Error: Cell cycle file missing: {CELL_CYCLE_FILE}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    print(f"[CopyVAE] Input: {h5ad_file}")
    print(f"[CopyVAE] Output: {args.outdir}")

    cmd = [
        "python", "-m", "copyvae.pipeline",
        h5ad_file,
        CELL_CYCLE_FILE,
        "--gene_map_csv", GENE_MAP_CSV,
    ]

    # Force CPU: hide the GPU so TF never grabs GPU memory / OOMs
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = ""
    try:
        subprocess.run(cmd, check=True, env=env)
        print("[CopyVAE] Pipeline finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"[CopyVAE] Failed with exit code {e.returncode}")
        sys.exit(1)

if __name__ == "__main__":
    main()