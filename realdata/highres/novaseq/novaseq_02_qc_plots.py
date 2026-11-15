#!/usr/bin/env python
"""Deliverable 1b — QC panels for the NovaSeq X Illumina spatial samples.

Reads the cell x gene h5ads from novaseq_01_prep_cells.py and writes, per sample,
nCount_spatial / nFeature_spatial distributions and spatial maps, plus a summary table.

Run:
  XDG_CACHE_HOME="$TMPDIR/fccache" \
  <CONDA_ROOT>/envs/bin2cell_env/bin/python \
    high_res_ST/novaseq/code/novaseq_02_qc_plots.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd

from novaseq_samples import SAMPLES, PROJ, OUTDIR

CELLS = os.path.join(PROJ, "processed_data", "cells")
MIN_COUNT, MIN_FEATURE = 50, 25   # candidate inferCNV inclusion thresholds


def save_both(fig, base):
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print("DONE ->", os.path.basename(base) + ".{png,pdf}", flush=True)


def spatial_panel(ax, df, col, title):
    v = df[col].values.astype(float)
    o = np.argsort(v)
    sc = ax.scatter(df.x.values[o], df.y.values[o], c=v[o], s=0.15, lw=0,
                    cmap="viridis", vmin=np.quantile(v, 0.01), vmax=np.quantile(v, 0.99))
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (um)", fontsize=8); ax.set_ylabel("y (um)", fontsize=8)
    ax.tick_params(labelsize=7)
    plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.02).ax.tick_params(labelsize=7)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    obs = {}
    for s in SAMPLES:
        p = os.path.join(CELLS, f"{s['name']}_cells.h5ad")
        if not os.path.exists(p):
            print("MISSING, skipping:", p, flush=True)
            continue
        obs[s["name"]] = ad.read_h5ad(p, backed="r").obs.copy()
        print(f"loaded {s['name']}: {len(obs[s['name']]):,} cells", flush=True)
    if not obs:
        raise SystemExit("no h5ads found - run novaseq_01_prep_cells.py first")

    # ---- summary table ----
    rows = []
    for name, df in obs.items():
        keep = (df.nCount_spatial >= MIN_COUNT) & (df.nFeature_spatial >= MIN_FEATURE)
        rows.append(dict(
            sample=name, n_cells=len(df),
            median_nCount=float(df.nCount_spatial.median()),
            mean_nCount=float(df.nCount_spatial.mean()),
            median_nFeature=float(df.nFeature_spatial.median()),
            mean_nFeature=float(df.nFeature_spatial.mean()),
            median_area_um2=float(df.area_um2.median()),
            median_barcodes_per_cell=float(df.n_barcodes.median()),
            total_counts_in_cells=int(df.nCount_spatial.sum()),
            n_cells_pass=int(keep.sum()),
            frac_cells_pass=round(float(keep.mean()), 3)))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(OUTDIR, "qc01_summary.csv"), index=False)
    print(summ.to_string(index=False), flush=True)

    # ---- distributions ----
    fig, axes = plt.subplots(len(obs), 3, figsize=(13, 3.6 * len(obs)))
    axes = np.atleast_2d(axes)
    for i, (name, df) in enumerate(obs.items()):
        for j, (col, thr, lab) in enumerate([
                ("nCount_spatial", MIN_COUNT, "counts per cell"),
                ("nFeature_spatial", MIN_FEATURE, "genes per cell")]):
            v = df[col].values
            axes[i, j].hist(v, bins=120, range=(0, np.quantile(v, 0.995)), color="#4C72B0")
            axes[i, j].axvline(thr, color="#C44E52", ls="--", lw=1)
            axes[i, j].set_xlabel(f"{col}  (x-axis clipped at p99.5)", fontsize=8)
            axes[i, j].set_ylabel("cells", fontsize=8)
            axes[i, j].set_title(f"{name} — {lab}\nmedian {df[col].median():.0f}", fontsize=9)
            axes[i, j].tick_params(labelsize=7)
        axes[i, 2].scatter(df.area_um2.values, df.nCount_spatial.values,
                           s=0.5, lw=0, alpha=0.15, color="#55A868")
        axes[i, 2].set_ylim(0, np.quantile(df.nCount_spatial.values, 0.995))
        axes[i, 2].set_xlim(0, np.quantile(df.area_um2.values, 0.995))
        axes[i, 2].set_xlabel("cell area (um2)", fontsize=8)
        axes[i, 2].set_ylabel("nCount_spatial", fontsize=8)
        axes[i, 2].set_title(f"{name} — counts vs area", fontsize=9)
        axes[i, 2].tick_params(labelsize=7)
    fig.tight_layout()
    save_both(fig, os.path.join(OUTDIR, "qc01_distributions"))

    # ---- violins across samples ----
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    names = list(obs)
    for ax, col in zip(axes, ["nCount_spatial", "nFeature_spatial"]):
        vals = [obs[n][col].values for n in names]
        ax.violinplot(vals, showmedians=True, widths=0.85)
        ax.set_ylim(0, max(np.quantile(v, 0.99) for v in vals))
        ax.set_xticks(range(1, len(names) + 1))
        ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel(f"{col}  (y-axis clipped at p99)", fontsize=9)
        ax.tick_params(labelsize=7)
    fig.tight_layout()
    save_both(fig, os.path.join(OUTDIR, "qc01_violin"))

    # ---- spatial maps ----
    fig, axes = plt.subplots(len(obs), 2, figsize=(11, 5.5 * len(obs)))
    axes = np.atleast_2d(axes)
    for i, (name, df) in enumerate(obs.items()):
        spatial_panel(axes[i, 0], df, "nCount_spatial", f"{name} — nCount_spatial")
        spatial_panel(axes[i, 1], df, "nFeature_spatial", f"{name} — nFeature_spatial")
    fig.tight_layout()
    save_both(fig, os.path.join(OUTDIR, "qc01_spatial"))


if __name__ == "__main__":
    main()
