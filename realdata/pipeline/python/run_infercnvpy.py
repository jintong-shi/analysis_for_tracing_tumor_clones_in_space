#!/usr/bin/env python3
"""
run_infercnvpy.py — Reproducible CLI wrapper for infercnvpy on simulated (or real) ST data

Features
- Accepts either an input .h5ad OR a gene-by-cell (spot) expression Parquet + annotation Parquet.
- Flexible reference definition: plain text list of reference spot IDs OR annotation + a reference group label.
- Merges a gene order file (gene, chromosome, start, end), harmonizes chromosome naming, and excludes chrM by default.
- Normalizes (CPM to target_sum) + log1p.
- Runs infercnvpy, saves results to .h5ad, CNV (gene-level) matrix to Parquet, and an optional chromosome heatmap.

Example usages
1) H5AD + reference ID list
   python run_infercnvpy.py \
     --input_h5ad /path/to/obj.h5ad \
     --ref_spots_file /path/to/ref_spots.txt \
     --gene_order_file /path/to/gencode_v38_gene_pos.txt \
     --out_dir ./infercnvpy_out --window_size 100 --n_jobs 8

2) Parquet matrices (sim data layout shown in your R script)
   python run_infercnvpy.py \
     --data_dir /path/to/sim_run \
     --gene_order_file /path/to/gencode_v38_gene_pos.txt \
     --ref_annotation_file /path/to/anno.parquet \
     --ref_group Benign \
     --out_dir ./infercnvpy_out

3) Batch over many simulation runs (bash)
   for d in /data/sim/*; do \
     python run_infercnvpy.py --data_dir "$d" \
       --gene_order_file /path/gencode_v38_gene_pos.txt \
       --ref_group Benign --out_dir "$d/infercnvpy_out"; \
   done
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

# ---------------------------
# Utilities
# ---------------------------

def log(msg: str):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] {msg}")


def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)

# ---------- 辅助函数（放在文件顶部 import 之后，或此处上方） ----------
def _get_cnv_matrix_cellsxgenes(adata: sc.AnnData) -> np.ndarray:
    """返回 cells x genes 的 CNV 矩阵；若不存在则抛异常。"""
    if "gene_values_cnv" in adata.layers and adata.layers["gene_values_cnv"] is not None:
        M = np.asarray(adata.layers["gene_values_cnv"])
        if M.shape == (adata.n_obs, adata.n_vars):
            return M
        if M.shape == (adata.n_vars, adata.n_obs):
            return M.T
    if "X_cnv" in adata.obsm and adata.obsm["X_cnv"] is not None:
        M = np.asarray(adata.obsm["X_cnv"])
        if M.shape == (adata.n_obs, adata.n_vars):
            return M
    raise ValueError("CNV matrix not found in layers['gene_values_cnv'] or obsm['X_cnv'].")

def _compute_gene_order_by_coord(adata: sc.AnnData) -> list[str]:
    """按 chr/start/end 生成基因顺序（用于热图列顺序）。"""
    need = {"chromosome","start","end"}
    if not need.issubset(set(adata.var.columns)):
        raise ValueError("adata.var 需要含有 columns: 'chromosome','start','end'")
    var = adata.var.copy()
    chr_clean = var["chromosome"].astype(str).str.replace("^chr","",regex=True).str.upper()
    def _rank(x):
        if x == "X": return 23
        if x == "Y": return 24
        if x in ("M","MT"): return 25
        try: return int(x)
        except: return 99
    var["__chr_rank"] = chr_clean.map(_rank)
    return var.sort_values(["__chr_rank","start","end"]).index.to_list()

def _compute_cell_order(adata: sc.AnnData, groupby: str|None) -> list[str]:
    """基于 CNV（cells x genes）做层次聚类得到行顺序；若提供 groupby，则组内聚类，按分类顺序拼接。"""
    X = _get_cnv_matrix_cellsxgenes(adata)
    cells = list(adata.obs_names)

    def _cluster_order(idx: np.ndarray) -> list[int]:
        Xi = X[idx, :]
        try:
            from scipy.spatial.distance import pdist
            from scipy.cluster.hierarchy import linkage, leaves_list
            Z = linkage(pdist(Xi, metric="correlation"), method="average")
            leaf_local = leaves_list(Z)
            return list(idx[leaf_local])
        except Exception:
            # 回退：按 |CNV| 强度排序
            mag = np.nanmean(np.abs(Xi), axis=1)
            return list(idx[np.argsort(-mag)])

    if groupby and groupby in adata.obs:
        cats = pd.Categorical(adata.obs[groupby])
        out_idx = []
        for cat in cats.categories:
            idx = np.where(cats == cat)[0]
            if len(idx) == 0: 
                continue
            out_idx.extend(_cluster_order(idx))
        return [cells[i] for i in out_idx]
    else:
        idx_all = np.arange(adata.n_obs)
        return [cells[i] for i in _cluster_order(idx_all)]

def _pick_fig(axthing):
    """把 infercnvpy 的返回对象转换为 matplotlib Figure。"""
    if axthing is None:
        return plt.gcf()
    if hasattr(axthing, "figure"):
        return axthing.figure
    if isinstance(axthing, (list, tuple, np.ndarray)) and len(axthing) > 0:
        first = axthing.flat[0] if hasattr(axthing, "flat") else axthing[0]
        if hasattr(first, "figure"):
            return first.figure
    return plt.gcf()

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

    # Duplicate var_names break attach_gene_order's reindexing/merge.
    if not adata.var_names.is_unique:
        n_dup = int(adata.var_names.duplicated().sum())
        log(f"WARNING: {n_dup} duplicate gene names; calling var_names_make_unique().")
        adata.var_names_make_unique()

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

    # ---------- groupby 逻辑（优先用 ref；否则用 infercnvpy 的 PCA/邻接/Leiden 生成 cnv_leiden） ----------
    groupby_key = None
    if ref_cat:
        groupby_key = ref_key
        log(f"Will plot grouped by reference: groupby={groupby_key}")
    else:
        # 若未有 CNV 邻接，就在 CNV 表征上显式跑 PCA → neighbors
        # 有些版本的 infercnvpy 接口支持 use_rep/key_added，有些不支持；做兼容处理
        cnv.tl.pca(adata)
        # neighbors on CNV representation
        # 目标是在 adata.uns["cnv_neighbors"] / adata.obsp["cnv_connectivities"] 等下生成 CNV 图
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata)
        groupby_key = "cnv_leiden"
        log(f"Computed cnv_leiden via infercnvpy; counts: {adata.obs['cnv_leiden'].value_counts().to_dict()}")

    # ---------- 计算并保存绘图顺序（基因/spot） ----------
    gene_order = None
    cell_order = None
    try:
        gene_order = _compute_gene_order_by_coord(adata)
        cell_order = _compute_cell_order(adata, groupby=groupby_key)
        gene_order_fp = os.path.join(args.out_dir, f"{args.out_prefix}_gene_order.txt")
        cell_order_fp = os.path.join(args.out_dir, f"{args.out_prefix}_cell_order.txt")
        pd.Series(gene_order, name="gene").to_csv(gene_order_fp, index=False)
        pd.Series(cell_order, name="cell").to_csv(cell_order_fp, index=False)
        log(f"Saved gene order → {gene_order_fp}")
        log(f"Saved cell order → {cell_order_fp}")
    except Exception as e:
        warnings.warn(f"Failed to compute/save plot orders: {e}")

    # ---------- 保存基因层面的 CNV 矩阵（genes x cells） ----------
    try:
        X = _get_cnv_matrix_cellsxgenes(adata)
        cnv_df = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
        parquet_fp = os.path.join(args.out_dir, f"{args.out_prefix}_gene_cnv.parquet")
        cnv_df.reset_index(names="gene").to_parquet(parquet_fp, compression="snappy")
        log(f"Saved gene-level CNV matrix → {parquet_fp}")
    except Exception as e:
        warnings.warn(f"Failed to extract CNV matrix: {e}")

    # ---------- 画热图（使用我们保存的顺序；无 ref 用 cnv_leiden 分组） ----------
    try:
        adata_plot = adata
        if cell_order is not None:
            keep_cells = [c for c in cell_order if c in adata.obs_names]
            adata_plot = adata_plot[keep_cells, :]
        if gene_order is not None:
            keep_genes = [g for g in gene_order if g in adata.var_names]
            adata_plot = adata_plot[:, keep_genes]

        ax_or_arr = cnv.pl.chromosome_heatmap(
            adata_plot,
            groupby=groupby_key if groupby_key else None,
            dendrogram=False  # 顺序由我们控制，便于复现
        )
        fig = _pick_fig(ax_or_arr)
        heatmap_fp = os.path.join(args.out_dir, f"{args.out_prefix}_heatmap.png")
        fig.savefig(heatmap_fp, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log(f"Saved heatmap → {heatmap_fp}")
    except Exception as e:
        warnings.warn(f"Failed to plot heatmap: {e}")

    # ---------- 写回 h5ad，并把 cnv_leiden 另存一份表 ----------
    out_h5ad = os.path.join(args.out_dir, f"{args.out_prefix}.h5ad")
    adata.write(out_h5ad)
    log(f"Saved inferCNV result → {out_h5ad}")
    if "cnv_leiden" in adata.obs.columns:
        leiden_csv = os.path.join(args.out_dir, f"{args.out_prefix}_cnv_leiden.csv")
        adata.obs[["cnv_leiden"]].to_csv(leiden_csv)
        log(f"Saved cnv_leiden table → {leiden_csv}")
    
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
