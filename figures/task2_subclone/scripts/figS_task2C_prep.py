#!/usr/bin/env python
"""One-time provenance step for figS_task2C (Task-2 subclone CNV heatmaps).

Slices the raw per-tool spots x genes matrices (~740 MB, live only in
`syntheticdata_cnv_profile/`) down to exactly the plotted spots x genes and
writes compact inputs into `input_data/synth_cnv_heatmap/`, so that the figure
script (`figS_task2C_cnv_heatmap.py`) replots off this folder alone.

To guarantee the compact figure is identical to the frozen original, the raw
loaders and the panel table are imported verbatim from the original script,
vendored alongside this one as `vendor/subclone_heatmap.py` (plus the
`vendor/figure_style.py` it imports); this step just applies them and saves the
resulting arrays.

Run on a machine that has the source tree (st-a100 / st-nexus):
  calicost_env/bin/python scripts/figS_task2C_prep.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ORIG = HERE / "vendor" / "subclone_heatmap.py"
OUT_DIR = HERE.parent / "input_data" / "synth_cnv_heatmap"


def _load_original():
    spec = importlib.util.spec_from_file_location("subclone_heatmap_orig", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    orig = _load_original()
    rep_dir = Path(orig.DEFAULT_REP)
    syn_dir = Path(orig.DEFAULT_SYN)
    rep_name = rep_dir.name
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Plotted gene order and spot order (verbatim from the original main()).
    gene_pos = orig.load_gene_order(syn_dir)
    genes = gene_pos["gene"].tolist()
    gt_hard = pd.read_csv(rep_dir / "subclone_groundtruth_spot_hard.csv")
    gt_hard = gt_hard.sort_values(["subclone", "spot_id"], kind="stable")
    spots = gt_hard["spot_id"].tolist()
    gt_lab = gt_hard["subclone"].tolist()

    gene_pos.to_parquet(OUT_DIR / "gene_order.parquet")
    pd.DataFrame({"spot_id": spots, "subclone": gt_lab}).to_csv(
        OUT_DIR / "spots.csv", index=False)

    arrays: dict[str, np.ndarray] = {}
    pred_cols: dict[str, pd.Series] = {}
    manifest = []

    # Ground truth panel (no predicted-subclone bar).
    gt_cnv = orig.load_gt_cnv(syn_dir, genes, spots).astype(np.float32)
    arrays["gt"] = gt_cnv
    manifest.append({"name": "Ground truth", "key": "gt", "has_pred": False,
                     "present": True})

    for name, kind, cnv_dir, sub_dir in orig.PANELS:
        key = cnv_dir
        if kind == "state":
            cnv = orig.load_infercnv_state(rep_dir, cnv_dir, genes, spots)
        else:
            cnv = orig.load_tool_cnv(rep_dir, cnv_dir, genes, spots)
        present = cnv is not None
        if present:
            arrays[key] = cnv.astype(np.float32)
        pred = orig.load_pred_labels(rep_dir, sub_dir, spots)
        pred_cols[name] = pred.reset_index(drop=True)
        manifest.append({"name": name, "key": key, "has_pred": True,
                         "present": present})

    np.savez_compressed(OUT_DIR / "cnv_arrays.npz", **arrays)
    pd.DataFrame(pred_cols).to_csv(OUT_DIR / "pred_labels.csv", index=False)
    (OUT_DIR / "manifest.json").write_text(json.dumps(
        {"rep_name": rep_name, "n_spots": len(spots), "n_genes": len(genes),
         "panels": manifest}, indent=2))

    total = sum(f.stat().st_size for f in OUT_DIR.iterdir() if f.is_file())
    print(f"[prep] rep={rep_name}  spots={len(spots)}  genes={len(genes)}")
    print(f"[prep] wrote {OUT_DIR}  ({total/1e6:.1f} MB)")
    for f in sorted(OUT_DIR.iterdir()):
        print(f"        {f.name:20s} {f.stat().st_size/1e6:8.2f} MB")


if __name__ == "__main__":
    main()
