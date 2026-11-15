#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, argparse
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc

def read_keep(path):
    with open(path) as f:
        lst = [x.strip() for x in f if x.strip()]
    return lst if lst else None

def align_keep_to_obs(obs_names, keep):
    """尽量对齐 keep 与数据条形码（处理 -1 后缀差异）；keep=None 则返回 None 表示保留全部。"""
    if not keep:
        return None
    obs_set = set(obs_names)
    direct = [k for k in keep if k in obs_set]
    if direct:
        return direct
    # 尝试加 -1
    alt = [k+"-1" for k in keep if (k+"-1") in obs_set]
    if alt:
        return alt
    # 尝试去 -1
    alt2 = [k[:-2] if k.endswith("-1") else k for k in keep]
    alt2 = [k for k in alt2 if k in obs_set]
    if alt2:
        return alt2
    raise SystemExit("[copyVAE-prepare] keep_barcodes 与数据条形码无交集（检查 -1 后缀）。")

def ensure_copyvae_columns(adata: ad.AnnData):
    """为 copyVAE 补齐必要列；以 var_names 为 gene 名称来源。"""
    g = adata.var_names.astype(str)
    if adata.var is None or adata.var.shape[1] == 0:
        adata.var = pd.DataFrame(index=g)
    # 优先使用 read_visium/read_10x_* 提供的 gene 名称列
    # 常见：var['gene_symbols'] 或 var['gene_names']
    if 'gene_symbols' in adata.var.columns:
        base = adata.var['gene_symbols'].astype(str)
    elif 'gene_names' in adata.var.columns:
        base = adata.var['gene_names'].astype(str)
    else:
        base = g
    adata.var['gene_name']       = base
    adata.var['gene_symbols']    = base
    adata.var['symbol']          = base
    adata.var['GeneName']        = base
    adata.var['name']            = base
    adata.var['gene_short_name'] = base

    # gene_id 相关（10x 通常有 gene_ids 列）
    if 'gene_ids' in adata.var.columns:
        gid = adata.var['gene_ids'].astype(str)
    elif 'gene_id' in adata.var.columns:
        gid = adata.var['gene_id'].astype(str)
    else:
        gid = g
    adata.var['gene_ids'] = gid
    adata.var['gene_id']  = gid
    adata.var['GeneID']   = gid

    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--visium_folder", required=True, help="通常是 SpaceRanger 的 outs 目录")
    ap.add_argument("--keep_barcodes_file", required=True)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    keep = read_keep(args.keep_barcodes_file)

    # 1) 优先尝试 read_visium（会自动寻找 outs 下的 h5）
    adata = None
    try:
        # load_images=False 可以避开图片依赖
        adata = sc.read_visium(args.visium_folder, load_images=False)
        print("[copyVAE-prepare] read_visium 成功")
    except Exception as e:
        print(f"[copyVAE-prepare] read_visium 失败，将回退到 read_10x_mtx: {e}")

    # 2) 回退到 read_10x_mtx（filtered_feature_bc_matrix 目录）
    if adata is None:
        mtx_dir = os.path.join(args.visium_folder, "filtered_feature_bc_matrix")
        if not os.path.isdir(mtx_dir):
            # 有些人传的就是 mtx 目录本身
            if os.path.isfile(os.path.join(args.visium_folder, "matrix.mtx")) or \
               os.path.isfile(os.path.join(args.visium_folder, "matrix.mtx.gz")):
                mtx_dir = args.visium_folder
            else:
                raise SystemExit(f"[copyVAE-prepare] 未找到 10x mtx 目录: {mtx_dir}")
        # 用 gene_symbols 更符合 copyVAE 的期望（若不存在会自动 fallback）
        adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
        if "gene_name" not in adata.var.columns:
                adata.var["gene_name"] = adata.var_names.astype(str)
        print(f"[copyVAE-prepare] read_10x_mtx 成功: {mtx_dir}")

    # 3) 过滤 keep barcodes（如提供）
    aligned = align_keep_to_obs(adata.obs_names.tolist(), keep)
    if aligned is not None:
        adata = adata[aligned, :].copy()

    # 4) 补齐 copyVAE 需要的列
    ensure_copyvae_columns(adata)

    # 5) 基本校验并写出
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise SystemExit(f"[copyVAE-prepare] 结果为空：cells={adata.n_obs}, genes={adata.n_vars}")

    # ---- 统一并补齐基因名列，确保 copyVAE 能识别 ----
    # 1) 先保证 var_names 唯一
    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()

    # 2) 如果存在常见列，优先使用；否则就用行名作为基因名
    preferred = None
    for cand in ["gene_name", "gene_symbols", "GeneName", "symbol", "gene", "GeneSymbol"]:
        if cand in adata.var.columns:
            preferred = cand
            break
    if preferred is None:
        preferred = "__from_index__"
        adata.var[preferred] = adata.var_names.astype(str)

    # 3) 写全所有可能别名，全部指向同一份干净字符串（避免类别/空值）
    vals = adata.var[preferred].astype(str).fillna("")
    for alias in ["gene_name", "gene_symbols", "GeneName", "GeneSymbol", "gene", "symbol"]:
        adata.var[alias] = vals

    # 4) 防御性地去掉首尾空白
    adata.var["gene_name"] = adata.var["gene_name"].str.strip()
    adata.var["symbol"]    = adata.var["symbol"].str.strip()

    # 5) 再保证一次唯一（某些 10x 符号会重复）
    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()

    out_h5ad = os.path.join(args.output_dir, "st_expression.h5ad")
    adata.write_h5ad(out_h5ad)
    with open(os.path.join(args.output_dir, "barcodes_order.txt"), "w") as f:
        for bc in adata.obs_names.tolist():
            f.write(bc + "\n")

    print(f"[copyVAE-prepare] Wrote {out_h5ad}  cells={adata.n_obs}  genes={adata.n_vars}")
    print(adata)

if __name__ == "__main__":
    sys.exit(main())