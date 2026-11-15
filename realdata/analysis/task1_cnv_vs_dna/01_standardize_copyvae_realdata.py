#!/usr/bin/env python3
"""Standardize copyVAE outputs into spots x genes log2 CSV.

Split out from 01_standardize_cnv_output_realdata.R because RcppCNPy::npyLoad
segfaults at the C level on some copy.npy files and aborts the entire R run.
Loading .npy via numpy here is segfault-safe, and a failure only skips that
one tool directory.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

EPS = 1e-12
DEFAULT_DATASETS_DIR = "<DATA_ROOT>/realdata/datasets"
RUN_DIR_CANDIDATES = ("run_CNV", "run_links")
TRUTH_SUBDIRS = ("wes_cnv", "wgs_cnv")
COPYVAE_DIRNAME_RE = re.compile(r"^copyvae$", re.IGNORECASE)


def log_msg(msg: str) -> None:
    print(msg, flush=True)


def normalize_key(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def normalize_keys(xs):
    return {normalize_key(x) for x in xs}


def normalize_spot_ids(ids):
    return [str(s).strip().replace(".", "-") for s in ids]


def parse_csv_arg(x):
    if x is None or x == "":
        return None
    out = [s.strip() for s in x.split(",")]
    out = [s for s in out if s]
    return out or None


def write_standard_csv(df: pd.DataFrame, out_file: Path) -> None:
    # Match R write.csv(mat, quote=FALSE): leading empty header cell, no quoting.
    df.to_csv(out_file, index=True, index_label="", quoting=3)  # 3 = csv.QUOTE_NONE


def load_copy_npy(path: Path):
    """Load a numpy .npy file safely. Returns 2D float ndarray or None."""
    if not path.is_file():
        return None
    try:
        arr = np.load(path, allow_pickle=False)
    except Exception as e:
        log_msg(f"[copyVAE] Skip: failed to load {path}: {e}")
        return None
    if arr.ndim != 2 or arr.size == 0:
        log_msg(f"[copyVAE] Skip: copy.npy must be a non-empty 2D matrix in {path.parent}")
        return None
    return arr.astype(np.float64, copy=False)


def run_copyvae(tool_dir: Path, eps: float, force: bool) -> bool:
    out_file = tool_dir / "processed_log2_spots_x_genes.csv"
    if out_file.exists() and not force:
        log_msg(f"[copyVAE] Skip existing: {out_file}")
        return True

    f_copy = tool_dir / "copy.npy"
    f_map = tool_dir / "bin_to_gene_map.csv"
    f_barc = tool_dir / "barcodes_order.txt"

    if not f_copy.is_file() or not f_map.is_file():
        log_msg(f"[copyVAE] Skip: copy.npy or bin_to_gene_map.csv missing in {tool_dir}")
        return False

    copy_mat = load_copy_npy(f_copy)
    if copy_mat is None:
        return False

    # copy.npy is spots x bins; bins_x_spots with baseline 2.
    bins_x_spots = (copy_mat / 2.0).T  # shape: (n_bins, n_spots)

    try:
        mapping = pd.read_csv(f_map)
    except Exception as e:
        log_msg(f"[copyVAE] Skip: failed to read {f_map}: {e}")
        return False
    if not {"gene", "bin_index"}.issubset(mapping.columns):
        log_msg(f"[copyVAE] Skip: bin_to_gene_map.csv must contain gene/bin_index in {tool_dir}")
        return False

    bin_idx = pd.to_numeric(mapping["bin_index"], errors="coerce")
    n_bins = bins_x_spots.shape[0]
    ok = bin_idx.notna() & (bin_idx >= 0) & (bin_idx < n_bins)
    n_drop = int((~ok).sum())
    if n_drop > 0:
        log_msg(f"[copyVAE] Warn: dropped {n_drop} out-of-range bin_index rows in {f_map}")
    if not ok.any():
        log_msg(f"[copyVAE] Skip: no valid bin_index rows in {f_map}")
        return False
    mapping = mapping.loc[ok].reset_index(drop=True)
    ridx = bin_idx.loc[ok].astype(int).to_numpy()

    genes = mapping["gene"].astype(str).str.strip()
    keep_gene = genes.ne("") & genes.ne("nan")
    if not keep_gene.any():
        log_msg(f"[copyVAE] Skip: no non-empty gene names in {f_map}")
        return False
    genes = genes.loc[keep_gene].to_numpy()
    ridx = ridx[keep_gene.to_numpy()]

    gene_x_spots = bins_x_spots[ridx, :]
    gene_x_spots_df = pd.DataFrame(gene_x_spots, index=genes)

    if gene_x_spots_df.index.duplicated().any():
        gene_x_spots_df = gene_x_spots_df.groupby(level=0, sort=False).mean()

    # log2(x + eps) — matches R's `log2(gene_x_spots + eps)`.
    with np.errstate(invalid="ignore", divide="ignore"):
        log2_vals = np.log2(gene_x_spots_df.to_numpy() + eps)
    final_mat = pd.DataFrame(
        log2_vals.T, index=gene_x_spots_df.columns, columns=gene_x_spots_df.index
    )

    if f_barc.is_file():
        with open(f_barc, encoding="utf-8") as fh:
            barcodes = [b.strip() for b in fh if b.strip()]
        if len(barcodes) == final_mat.shape[0]:
            final_mat.index = normalize_spot_ids(barcodes)
        else:
            log_msg(
                f"[copyVAE] Warn: barcodes length ({len(barcodes)}) != nspots "
                f"({final_mat.shape[0]}). Keep existing rownames."
            )
            final_mat.index = normalize_spot_ids(final_mat.index)
    else:
        final_mat.index = normalize_spot_ids(final_mat.index)

    bad = [r for r in final_mat.index if r is None or r == ""]
    if bad:
        log_msg(f"[copyVAE] Skip: missing spot barcodes after processing in {tool_dir}")
        return False

    write_standard_csv(final_mat, out_file)
    log_msg(f"[copyVAE] Saved: {out_file}")
    return True


def find_run_dir(dataset_dir: Path):
    for cand in RUN_DIR_CANDIDATES:
        p = dataset_dir / cand
        if p.is_dir():
            return p
    return None


def discover_dataset_roots(datasets_dir: Path, datasets_filter=None):
    out = {}
    filter_keys = normalize_keys(datasets_filter) if datasets_filter else None
    for d in sorted(p for p in datasets_dir.iterdir() if p.is_dir()):
        if filter_keys is not None and normalize_key(d.name) not in filter_keys:
            continue
        truth_dir = d / "truth_links"
        has_truth = any((truth_dir / s).is_dir() for s in TRUTH_SUBDIRS)
        if not has_truth:
            continue
        run_dir = find_run_dir(d)
        if run_dir is None:
            log_msg(f"[Init] Skip dataset {d.name}: no {'/'.join(RUN_DIR_CANDIDATES)} present")
            continue
        out[d.name] = run_dir
    return out


def process_root(root_dir: Path, dataset_name: str, samples, strategies, eps, force):
    sample_keys = normalize_keys(samples) if samples else None
    strat_keys = normalize_keys(strategies) if strategies else None

    sample_dirs = sorted(p for p in root_dir.iterdir() if p.is_dir())
    if sample_keys is not None:
        sample_dirs = [p for p in sample_dirs if normalize_key(p.name) in sample_keys]

    log_msg(f"\n[Dataset] {dataset_name} | root_dir = {root_dir} | samples = {len(sample_dirs)}")

    for sample_dir in sample_dirs:
        run_dirs = sorted(p for p in sample_dir.iterdir() if p.is_dir())
        if strat_keys is not None:
            run_dirs = [p for p in run_dirs if normalize_key(p.name) in strat_keys]
        log_msg(f"\n[Sample] {sample_dir.name} | runs = {len(run_dirs)}")

        for run_dir in run_dirs:
            copyvae_dirs = [
                p for p in run_dir.iterdir()
                if p.is_dir() and COPYVAE_DIRNAME_RE.match(p.name)
            ]
            log_msg(f"  [Run] {run_dir.name} | copyVAE dirs = {len(copyvae_dirs)}")
            for tool_dir in copyvae_dirs:
                log_msg(f"    [Tool] {tool_dir.name} (copyvae)")
                try:
                    run_copyvae(tool_dir, eps=eps, force=force)
                except Exception as e:
                    log_msg(f"    [Tool] ERROR in {tool_dir.name}: {e}")
                    traceback.print_exc()


def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="Standardize copyVAE outputs (Python replacement for the R copyVAE branch)."
    )
    ap.add_argument("--datasets_dir", default=DEFAULT_DATASETS_DIR)
    ap.add_argument("--root_dir", default=None,
                    help="Single dataset's run dir; overrides --datasets_dir/--datasets.")
    ap.add_argument("--datasets", default=None, help="Comma-separated dataset names.")
    ap.add_argument("--samples", default=None, help="Comma-separated sample names.")
    ap.add_argument("--strategies", default=None,
                    help="Comma-separated strategy/run names (e.g. AllSpots_no-ref).")
    ap.add_argument("--eps", type=float, default=EPS)
    ap.add_argument("--force", action="store_true")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not (np.isfinite(args.eps) and args.eps > 0):
        raise SystemExit("--eps must be > 0")

    datasets = parse_csv_arg(args.datasets)
    samples = parse_csv_arg(args.samples)
    strategies = parse_csv_arg(args.strategies)

    if args.root_dir is not None:
        root_dir = Path(args.root_dir)
        if not root_dir.is_dir():
            raise SystemExit(f"root_dir does not exist: {root_dir}")
        # Match R: dataset name = parent-of-parent (dataset dir) basename.
        ds_name = root_dir.parent.name
        roots = {ds_name: root_dir}
    else:
        datasets_dir = Path(args.datasets_dir)
        if not datasets_dir.is_dir():
            raise SystemExit(f"datasets_dir does not exist: {datasets_dir}")
        roots = discover_dataset_roots(datasets_dir, datasets_filter=datasets)
        if not roots:
            raise SystemExit(f"No datasets with WES/WGS truth found under: {datasets_dir}")

    log_msg(f"[Init] datasets = {len(roots)} ({', '.join(roots)})")
    if strategies:
        log_msg(f"[Init] strategies filter = {', '.join(strategies)}")
    log_msg(f"[Init] force = {str(args.force).upper()}")

    for ds_name, root in roots.items():
        process_root(root, ds_name, samples=samples, strategies=strategies,
                     eps=args.eps, force=args.force)

    log_msg("\n[Done] copyVAE standardization finished.")


if __name__ == "__main__":
    main()
