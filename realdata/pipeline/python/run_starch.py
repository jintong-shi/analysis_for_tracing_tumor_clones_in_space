#!/usr/bin/env python3
"""
Prepare STARCH input from Visium (Space Ranger) outputs.

Default: build a subset 10X directory containing
  barcodes.tsv, features.tsv, matrix.mtx, tissue_positions_list.csv
all in the SAME FOLDER, filtered by keep_barcodes.txt, then call STARCH.

Option: --as_csv to write Gene x Spot CSV with first row "rowxcol" headers.

Usage (10X dir mode):
  python run_starch.py \
    --visium_folder <spaceranger outs folder or its parent> \
    --output_dir <out> \
    --keep_barcodes_file keep_barcodes.txt \
    --n_clusters 3 --output_name STARCH
"""

import argparse, os, sys, subprocess, gzip, shutil
import pandas as pd
import numpy as np
from scipy.io import mmread, mmwrite
from scipy import sparse
from pathlib import Path

def first_existing(paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None

def require(path, msg):
    if not path or not os.path.exists(path):
        raise FileNotFoundError(msg + f": {path}")

def read_lines(path):
    with open(path, "r") as f:
        return [ln.strip() for ln in f if ln.strip()]

def expand_h5_to_mtx(h5_path: Path, out_dir: Path) -> Path:
    """Read a 10x .h5 (CellRanger v2 or v3+) and write an MTX-style triplet
    (matrix.mtx, barcodes.tsv, features.tsv) into out_dir."""
    import h5py
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _decode(arr):
        return [x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]

    with h5py.File(str(h5_path), "r") as f:
        if "matrix" in f:  # CellRanger v3+
            grp = f["matrix"]
            barcodes = _decode(grp["barcodes"][:])
            data = grp["data"][:]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]
            shape = tuple(int(x) for x in grp["shape"][:])
            feat = grp["features"]
            feat_ids = _decode(feat["id"][:])
            feat_names = _decode(feat["name"][:])
            feat_types = _decode(feat["feature_type"][:]) if "feature_type" in feat \
                else ["Gene Expression"] * len(feat_ids)
        else:  # CellRanger v2 (genome-named group)
            keys = [k for k in f.keys() if isinstance(f[k], h5py.Group)]
            if not keys:
                raise ValueError(f"Unrecognized 10x .h5 layout: {h5_path}")
            grp = f[keys[0]]
            barcodes = _decode(grp["barcodes"][:])
            data = grp["data"][:]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]
            shape = tuple(int(x) for x in grp["shape"][:])
            feat_ids = _decode(grp["genes"][:])
            feat_names = _decode(grp["gene_names"][:])
            feat_types = ["Gene Expression"] * len(feat_ids)

    M = sparse.csc_matrix((data, indices, indptr), shape=shape)  # genes x cells
    mmwrite(str(out_dir / "matrix.mtx"), M.tocoo())
    with open(out_dir / "barcodes.tsv", "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with open(out_dir / "features.tsv", "w") as fh:
        for i, n, t in zip(feat_ids, feat_names, feat_types):
            fh.write(f"{i}\t{n}\t{t}\n")
    return out_dir

def load_10x_handles(visium_folder, h5_stage_dir=None):
    """
    Return dict with:
      tenx_dir: filtered_feature_bc_matrix directory (or a directory that contains matrix.mtx*)
      barcodes, features, matrix: file paths
      tpos: tissue_positions_list.csv or tissue_positions.csv path
    Tries common layouts:
      <visium>/outs/filtered_feature_bc_matrix + <visium>/outs/spatial/...
      <visium>/filtered_feature_bc_matrix     + <visium>/spatial/...
      <visium> (already the matrix dir)       + sibling/../spatial/...
    Falls back to filtered_feature_bc_matrix.h5 if no MTX directory is present
    (requires h5_stage_dir to write the expanded triplet to).
    """
    vf = Path(visium_folder)

    # candidates for matrix dir
    candidates = [
        vf / "outs" / "filtered_feature_bc_matrix",
        vf / "filtered_feature_bc_matrix",
        vf,  # in case visium_folder is already the matrix dir
    ]
    tenx = None
    for cand in candidates:
        if (cand / "matrix.mtx").is_file() or (cand / "matrix.mtx.gz").is_file():
            tenx = cand
            break

    if not tenx:
        # Fallback: expand filtered_feature_bc_matrix.h5 if available
        h5_candidates = [
            vf / "outs" / "filtered_feature_bc_matrix.h5",
            vf / "filtered_feature_bc_matrix.h5",
        ]
        h5_path = next((p for p in h5_candidates if p.is_file()), None)
        if h5_path is None:
            raise FileNotFoundError(
                "Cannot find filtered_feature_bc_matrix (matrix.mtx[.gz]) "
                "or filtered_feature_bc_matrix.h5 under visium_folder."
            )
        if h5_stage_dir is None:
            raise ValueError("h5_stage_dir must be provided to expand .h5 input.")
        print(f"[STARCH] Expanding {h5_path} -> {h5_stage_dir}")
        tenx = expand_h5_to_mtx(h5_path, Path(h5_stage_dir))

    bar = first_existing([str(tenx / "barcodes.tsv.gz"), str(tenx / "barcodes.tsv")])
    fea = first_existing([
        str(tenx / "features.tsv.gz"),
        str(tenx / "features.tsv"),
        str(tenx / "genes.tsv.gz"),
        str(tenx / "genes.tsv"),
    ])
    mat = first_existing([str(tenx / "matrix.mtx.gz"), str(tenx / "matrix.mtx")])
    require(bar, "Missing barcodes.tsv(.gz)")
    require(fea, "Missing features.tsv(.gz)/genes.tsv(.gz)")
    require(mat, "Missing matrix.mtx(.gz)")

    # tissue positions (old & new names) — anchor on visium_folder, not tenx,
    # because tenx may live under output_dir when expanded from .h5
    spatial_candidates = [
        vf / "outs" / "spatial",
        vf / "spatial",
        tenx.parent / "spatial",
        tenx.parent.parent / "spatial",
    ]
    tpos = None
    for sp in spatial_candidates:
        if sp and sp.is_dir():
            tpos = first_existing([str(sp / "tissue_positions_list.csv"), str(sp / "tissue_positions.csv")])
            if tpos:
                break
    require(tpos, "Missing tissue positions (tissue_positions_list.csv or tissue_positions.csv)")

    return dict(tenx_dir=str(tenx), barcodes=bar, features=fea, matrix=mat, tpos=tpos)

def read_barcodes(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return [ln.strip() for ln in f]
    return [ln.strip() for ln in open(path)]

def read_features_df(path):
    # usually no header; keep as-is; do not drop columns
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", header=None, compression="gzip")
    return pd.read_csv(path, sep="\t", header=None)

def read_matrix(path):
    return mmread(path)  # genes x cells (sparse)

def load_tissue_positions(tpos_path):
    # old file: no header, 6 cols
    # new file: header with: barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
    base = os.path.basename(tpos_path)
    if base == "tissue_positions_list.csv":
        df = pd.read_csv(tpos_path, header=None)
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    else:
        df = pd.read_csv(tpos_path, header=0)
        col_map = {
            "barcode": "barcode",
            "in_tissue": "in_tissue",
            "array_row": "array_row",
            "array_col": "array_col",
            "pxl_row_in_fullres": "pxl_row",
            "pxl_col_in_fullres": "pxl_col",
            "pxl_row": "pxl_row",
            "pxl_col": "pxl_col",
        }
        df = df.rename(columns={c: col_map.get(c, c) for c in df.columns})
        need = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
        missing = [c for c in need if c not in df.columns]
        if missing:
            raise ValueError(f"tissue_positions.csv missing columns: {missing}")
        df = df[need]
    return df

def subset_10x_and_stage(tenx_info, keep_bcs, out_dir):
    """Write a subset 10X folder containing the 4 files in the SAME directory."""
    stage = os.path.join(out_dir, "starch_10x_subset")
    os.makedirs(stage, exist_ok=True)

    # read originals
    barcodes = read_barcodes(tenx_info["barcodes"])
    idx_map = {bc:i for i,bc in enumerate(barcodes)}

    # 当 keep_bcs 为 None => 使用全部
    if keep_bcs is None:
        keep_idx = list(range(len(barcodes)))
    else:
        keep_idx = [idx_map[bc] for bc in keep_bcs if bc in idx_map]
        if len(keep_idx) == 0:
            raise ValueError("No overlap between keep_barcodes and 10X barcodes after harmonization.")

    features_df = read_features_df(tenx_info["features"])
    M = read_matrix(tenx_info["matrix"]).tocsr()  # genes x cells

    # subset matrix by columns (cells/spots)
    M_sub = M[:, keep_idx]
    # write matrix.mtx
    mmwrite(os.path.join(stage, "matrix.mtx"), M_sub.tocoo())

    # write barcodes.tsv
    with open(os.path.join(stage, "barcodes.tsv"), "w") as f:
        for i in keep_idx:
            f.write(barcodes[i] + "\n")

    # write features.tsv
    features_out = os.path.join(stage, "features.tsv")
    features_df.to_csv(features_out, sep="\t", header=False, index=False)

    # write tissue_positions_list.csv (filtered)
    tpos_df = load_tissue_positions(tenx_info["tpos"])
    tpos_sub = tpos_df[tpos_df["barcode"].isin([barcodes[i] for i in keep_idx])].copy()
    tpos_out = os.path.join(stage, "tissue_positions_list.csv")
    tpos_sub.to_csv(tpos_out, header=False, index=False)

    return stage

def build_csv_from_subset(tenx_stage_dir):
    """Optional: construct Gene x Spot CSV from staged 10X dir (for --as_csv)."""
    mat = mmread(os.path.join(tenx_stage_dir, "matrix.mtx")).tocsr()  # genes x spots
    genes = pd.read_csv(os.path.join(tenx_stage_dir, "features.tsv"), sep="\t", header=False).iloc[:,1].astype(str).tolist()
    bcs   = read_barcodes(os.path.join(tenx_stage_dir, "barcodes.tsv"))
    tpos  = load_tissue_positions(os.path.join(tenx_stage_dir, "tissue_positions_list.csv"))
    pos_map = dict(zip(tpos["barcode"], [f"{r}x{c}" for r,c in zip(tpos["array_row"], tpos["array_col"])]))
    headers = [pos_map.get(bc, bc) for bc in bcs]  # fallback to bc if coords missing

    df = pd.DataFrame.sparse.from_spmatrix(mat.T, index=headers, columns=genes)
    df = df.T  # genes x spots; columns are headers
    df.insert(0, "gene", df.index)
    return df

def ensure_starch_repo(repo_dir: str | None, repo_url: str):
    """
    确保能找到 STARCH 的 run_STARCH.py。
    - 优先在 repo_dir/、repo_dir/STARCH/、repo_dir/STARCH_repo/ 下寻找。
    - 若不存在，则执行 clone（或 git pull）。
    返回 (repo_root_dir, run_script_path)。
    """
    def find_script(base: Path) -> Path | None:
        candidates = [
            base / "run_STARCH.py",
            base / "STARCH" / "run_STARCH.py",
            base / "STARCH_repo" / "run_STARCH.py",
        ]
        for p in candidates:
            if p.is_file():
                return p
        return None

    default_root = Path("./software/STARCH").resolve()
    root = Path(repo_dir).resolve() if repo_dir else default_root

    # 1) 已存在就先找脚本
    script = find_script(root)
    if script:
        return (script.parent, script)

    # 2) 如果 root 是个 git 仓库，先尝试 pull
    if root.is_dir() and (root / ".git").is_dir():
        try:
            subprocess.run(["git", "-C", str(root), "pull", "--ff-only"], check=True)
        except subprocess.CalledProcessError:
            pass
        script = find_script(root)
        if script:
            return (script.parent, script)

    # 3) 选择 clone 目标目录
    if not root.exists():
        root.parent.mkdir(parents=True, exist_ok=True)
        clone_target = root
    else:
        clone_target = root / "STARCH"
        clone_target.parent.mkdir(parents=True, exist_ok=True)
        # 清理空壳
        if clone_target.exists():
            try:
                # 若是空目录，删掉；否则保留
                next(clone_target.iterdir())
            except StopIteration:
                shutil.rmtree(clone_target, ignore_errors=True)

    # 4) clone 或 pull
    if clone_target.exists() and (clone_target / ".git").is_dir():
        subprocess.run(["git", "-C", str(clone_target), "pull", "--ff-only"], check=True)
    else:
        subprocess.run(["git", "clone", repo_url, str(clone_target)], check=True)

    # 5) 验证脚本
    script = find_script(clone_target) or (clone_target / "run_STARCH.py")
    if not script.is_file():
        raise FileNotFoundError(f"run_STARCH.py not found under {clone_target}")

    return (script.parent, script)

def harmonize_keep_with_data(barcodes_all_list, keep_list_or_none):
    """
    根据数据条形码与 keep 的后缀差异（是否带 '-1'）做自动对齐。
    - keep 为 None/空文件 => 返回 None（表示使用全部）
    - 有交集 => 返回按 keep 顺序过滤后的列表
    - 无交集：
        * 若 keep 全部不带 '-1' 且数据大多带 '-1'，给 keep 自动加 '-1'
        * 若 keep 大多带 '-1' 且数据不带，给 keep 自动去 '-1'
      仍无交集则报错
    """
    if keep_list_or_none is None or len(keep_list_or_none) == 0:
        return None  # 使用全部

    data_set = set(barcodes_all_list)
    keep = list(dict.fromkeys([k for k in keep_list_or_none if k]))  # 去重保序
    # 直接交集
    direct = [k for k in keep if k in data_set]
    if len(direct) > 0:
        return direct

    # 尝试加 '-1'
    data_has_dash = any(b.endswith("-1") for b in data_set)
    keep_has_dash = any(k.endswith("-1") for k in keep)
    if not keep_has_dash and data_has_dash:
        keep_alt = [k + "-1" for k in keep]
        alt = [k for k in keep_alt if k in data_set]
        if len(alt) > 0:
            return alt

    # 尝试去 '-1'
    if keep_has_dash and not data_has_dash:
        keep_alt = [k[:-2] if k.endswith("-1") else k for k in keep]
        alt = [k for k in keep_alt if k in data_set]
        if len(alt) > 0:
            return alt

    raise ValueError(
        f"keep_barcodes has no overlap with the data barcodes "
        f"(n_data={len(data_set)}, n_keep={len(keep)}). "
        f"Check '-1' suffix or leave keep file empty to use all spots."
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--visium_folder", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--keep_barcodes_file", required=True)
    ap.add_argument("--n_clusters", type=int, default=3)
    ap.add_argument("--output_name", default="STARCH")
    ap.add_argument("--starch_repo_dir", default="")
    ap.add_argument("--starch_repo_url", default="https://github.com/raphael-group/STARCH.git")
    ap.add_argument("--as_csv", action="store_true", help="Write Gene x Spot CSV instead of 10X dir")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # locate 10X + tissue positions (with .h5 fallback under output_dir)
    h5_stage = os.path.join(args.output_dir, "h5_to_mtx")
    tenx = load_10x_handles(args.visium_folder, h5_stage_dir=h5_stage)

    # read all barcodes from data (保持原始顺序)
    bar_all_list = read_barcodes(tenx["barcodes"])

    # keep subset
    keep_list = None
    if args.keep_barcodes_file and os.path.isfile(args.keep_barcodes_file):
        with open(args.keep_barcodes_file) as fh:
            keep_list = [ln.strip() for ln in fh if ln.strip()]
        if len(keep_list) == 0:
            keep_list = None  # 空文件 => 使用全部

    # 和数据条形码对齐（自动处理 '-1'）
    keep_aligned = harmonize_keep_with_data(bar_all_list, keep_list)

    # stage 10X subset（None => 全部）
    stage_dir = subset_10x_and_stage(tenx, keep_aligned, args.output_dir)

    # prepare STARCH input path
    starch_input_path = stage_dir
    if args.as_csv:
        df = build_csv_from_subset(stage_dir)
        csv_generated = os.path.join(args.output_dir, "gene_expression_matrix.csv")
        df.to_csv(csv_generated, index=False)
        starch_input_path = csv_generated

    # STARCH repo
    repo_root, run_py = ensure_starch_repo(args.starch_repo_dir, args.starch_repo_url)
    print(f"[STARCH] Using repo at: {repo_root}")
    print(f"[STARCH] Using runner: {run_py}")

    # run STARCH（在仓库根目录执行，保证能找到 hgTables_hg19.txt 等资源）
    cmd = [
        sys.executable,
        str(run_py),
        "-i", starch_input_path,
        "--output", args.output_name,
        "--n_clusters", str(args.n_clusters),
        "--outdir", args.output_dir
    ]
    subprocess.run(cmd, check=True, cwd=str(repo_root))

    open(os.path.join(args.output_dir, "STARCH.OK"), "w").close()

if __name__ == "__main__":
    main()