#!/usr/bin/env python3
"""STARCH input prep: Parquet expression matrix + metadata -> the CSV STARCH reads."""

import os
import argparse
import numpy as np
import pandas as pd
import sys

def ensure_unique_genes(df_expr):
    """If gene names (index) are duplicated, aggregate by sum."""
    if df_expr.index.has_duplicates:
        print(f"[WARN] Found {df_expr.index.duplicated().sum()} duplicated gene names; aggregating by sum.")
        df_expr = df_expr.groupby(level=0).sum()
    return df_expr

def load_expr(expr_path: str) -> pd.DataFrame:
    print(f"[INFO] Reading expression: {expr_path}")
    if not os.path.exists(expr_path):
        raise FileNotFoundError(f"Expression file not found: {expr_path}")
        
    df = pd.read_parquet(expr_path)
    if "index" not in df.columns:
        raise ValueError("Expression parquet must have a column named 'index' with gene names.")
    df = df.set_index("index")
    return df

def load_meta(meta_path: str) -> pd.DataFrame:
    print(f"[INFO] Reading metadata:   {meta_path}")
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    df = pd.read_parquet(meta_path)
    for col in ("index", "x", "y"):
        if col not in df.columns:
            raise ValueError(f"Metadata parquet must contain column '{col}'.")
    df = df.set_index("index")[["x", "y"]]
    return df

def main():
    ap = argparse.ArgumentParser(description="Convert Parquet expression+meta to STARCH CSV.")
    ap.add_argument("--data_dir",  required=True, help="Parent dir containing spatial_data/")
    ap.add_argument("--work_dir",  required=True, help="Where to save the output CSV")
    ap.add_argument("--input_dir", default=None,
                   help="Optional parquet file with 'cell_id' column for spot filtering.")
    args = ap.parse_args()

    data_dir = args.data_dir
    work_dir = args.work_dir
    anno_file = args.input_dir

    os.makedirs(work_dir, exist_ok=True)

    # Define paths to raw data
    # Assumes structure: {data_dir}/spatial_data/st-expression.parquet
    expr_path = os.path.join(data_dir, "spatial_data", "st-expression.parquet")
    meta_path = os.path.join(data_dir, "spatial_data", "st-meta.parquet")

    # 1) Load Data
    try:
        df_expr = load_expr(expr_path)        # genes x spots
        df_meta = load_meta(meta_path)        # spots x (x,y)
    except Exception as e:
        print(f"[ERROR] Failed to load input files: {e}")
        sys.exit(1)

    # 2) Optional Filter by Annotation
    if anno_file:
        print(f"[INFO] Filtering spots using annotation: {anno_file}")
        if not os.path.exists(anno_file):
            print(f"[ERROR] Annotation file not found: {anno_file}")
            sys.exit(1)
            
        df_anno = pd.read_parquet(anno_file)
        if "cell_id" not in df_anno.columns:
            raise ValueError("Annotation file must contain a 'cell_id' column.")
        
        keep_spots = (
            df_anno["cell_id"].astype(str).drop_duplicates().tolist()
        )
        # Keep order provided by annotation, drop those not in meta
        keep_spots = [s for s in keep_spots if s in df_meta.index]
        
        if not keep_spots:
            raise ValueError("After filtering, no overlapping spots between annotation and metadata.")
        
        df_meta = df_meta.reindex(keep_spots)
        print(f"[INFO] Retained {len(df_meta)} spots after filtering.")

    # 3) Align Expression Columns to Metadata Index
    common = [s for s in df_meta.index if s in df_expr.columns]
    if not common:
        raise ValueError("No overlapping spots between metadata and expression.")
    
    if len(common) != len(df_meta.index):
        missing = set(df_meta.index) - set(common)
        print(f"[WARN] {len(missing)} spots present in meta but missing in expression; they will be dropped.")
        df_meta = df_meta.loc[common]
        
    df_expr = df_expr.loc[:, df_meta.index]

    # 4) Create STARCH Coordinate Names "<x>x<y>"
    # Round to nearest int to avoid float naming issues
    coord_names = df_meta.apply(lambda r: f"{int(round(r['x']))}x{int(round(r['y']))}", axis=1)
    
    # Handle duplicates (e.g., spots very close to each other rounding to same int)
    if coord_names.duplicated().any():
        dup_n = coord_names.duplicated().sum()
        print(f"[WARN] {dup_n} duplicated coordinate names detected. Appending unique suffix.")
        seen = {}
        new_names = []
        for name in coord_names.tolist():
            if name in seen:
                seen[name] += 1
                new_names.append(f"{name}_{seen[name]}")
            else:
                seen[name] = 0
                new_names.append(name)
        coord_names = pd.Index(new_names)

    # 5) Ensure Numeric Data
    # Coerce errors to NaN then fill with 0
    df_expr = df_expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df_expr = ensure_unique_genes(df_expr).astype(np.float32)

    # 6) Assign Final Columns and Write
    df_expr.columns = coord_names
    
    out_csv = os.path.join(work_dir, "gene_expression_matrix.csv")
    
    print(f"[INFO] Writing CSV to {out_csv}...")
    # index_label="Gene" is required by STARCH
    df_expr.to_csv(out_csv, index_label="Gene")
    
    print(f"[SUCCESS] Saved STARCH-ready CSV.")
    print(f"          Shape: {df_expr.shape[0]} genes x {df_expr.shape[1]} spots")

if __name__ == "__main__":
    main()