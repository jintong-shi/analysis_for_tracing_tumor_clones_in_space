#!/usr/bin/env python
"""figS_task2C -- Task-2 subclone ground-truth vs tool CNV heatmaps (synthetic).

Standalone Task-2 supplementary figure. A ground-truth CNV heatmap on top, then
each tool's predicted CNV heatmap stacked below; all panels share the spot
ordering (by ground-truth subclone) and the genomic gene order on x, so the
ground truth can be compared against each tool at a glance. Left bars per row =
ground-truth subclone | predicted subclone. InferCNV CNV = HMM17 TumorFreeNormal
state (loss/neutral/gain -> -1/0/+1 on the same RdBu scale).

Self-contained: reads only `input_data/synth_cnv_heatmap/` (compact slices
produced by `figS_task2C_prep.py`) and writes to `supp_figures/`. The style
helpers used from the original `scripts/figures/figure_style.py` are inlined
below (only the pieces this figure needs; no figure_settings.tsv dependency).

Rendered with calicost_env/bin/python (matplotlib 3.10.1), matching the frozen
original in `paper/figures/fig3_syn_task2_subclone_ari/`.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Nimbus Sans", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

# --- style constants + helpers (inlined verbatim from figure_style.py) --------
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CLIP = 1.0
SUBCLONE_GT_CMAP = "Set1"
SUBCLONE_PRED_CMAP = "tab20"
CNV_CMAP = "RdBu_r"
NO_CALL_COLOR = 0.85


def _sort_clones(labels: Iterable[object]) -> list[str]:
    return sorted(
        {str(v) for v in labels if pd.notna(v)},
        key=lambda v: (int(m.group(1)) if (m := re.search(r"(\d+)", v)) else 1e9, v),
    )


def clone_palette(labels: Iterable[object], cmap_name: str) -> dict[str, object]:
    clones = _sort_clones(labels)
    cmap = plt.get_cmap(cmap_name)
    n_colors = getattr(cmap, "N", 20)
    return {clone: cmap(i % n_colors) for i, clone in enumerate(clones)}


def bar_rgb(labels: Iterable[object], palette: dict[str, object],
            missing: float = NO_CALL_COLOR) -> np.ndarray:
    labels = list(labels)
    rgb = np.full((len(labels), 1, 3), missing, dtype=float)
    for i, label in enumerate(labels):
        text = str(label)
        if pd.notna(label) and text in palette:
            rgb[i, 0, :] = palette[text][:3]
    return rgb


def chrom_boundaries(gene_pos: pd.DataFrame):
    chroms = gene_pos["chr"].to_numpy()
    bounds = [0] + [i for i in range(1, len(chroms)) if chroms[i] != chroms[i - 1]]
    bounds.append(len(chroms))
    mids = [((bounds[i] + bounds[i + 1]) / 2, chroms[bounds[i]])
            for i in range(len(bounds) - 1)]
    return bounds[1:-1], mids
# -----------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
IN_DIR = HERE.parent / "input_data" / "synth_cnv_heatmap"
OUT_DIR = HERE.parent / "supp_figures"


def main():
    manifest = json.loads((IN_DIR / "manifest.json").read_text())
    rep_name = manifest["rep_name"]
    gene_pos = pd.read_parquet(IN_DIR / "gene_order.parquet")
    spots_df = pd.read_csv(IN_DIR / "spots.csv")
    gt_lab = spots_df["subclone"].tolist()
    arrays = np.load(IN_DIR / "cnv_arrays.npz")
    pred_df = pd.read_csv(IN_DIR / "pred_labels.csv")

    panels = []
    pred_all = []
    for p in manifest["panels"]:
        cnv = arrays[p["key"]] if p["present"] else None
        if p["has_pred"]:
            pred = pred_df[p["name"]]
            pred_all.append(pred)
        else:
            pred = None
        panels.append((p["name"], cnv, gt_lab, pred))

    gt_palette = clone_palette(gt_lab, SUBCLONE_GT_CMAP)
    pred_palette = clone_palette(
        [x for p in pred_all for x in p], SUBCLONE_PRED_CMAP)

    sep_idx, chr_mids = chrom_boundaries(gene_pos)
    grp_runs = [i for i in range(1, len(gt_lab)) if gt_lab[i] != gt_lab[i - 1]]
    n_panels = len(panels)

    fig = plt.figure(figsize=(18, 2.5 * n_panels + 1.2))
    gs = GridSpec(n_panels, 3, figure=fig, width_ratios=[0.5, 0.5, 50],
                  hspace=0.32, wspace=0.03, left=0.07, right=0.91,
                  top=0.94, bottom=0.06)
    cmap = plt.get_cmap(CNV_CMAP).copy()
    cmap.set_bad("white")

    im = None
    for row, (name, cnv, gbar, pbar) in enumerate(panels):
        ax_g = fig.add_subplot(gs[row, 0])
        ax_p = fig.add_subplot(gs[row, 1])
        ax_h = fig.add_subplot(gs[row, 2])
        is_bottom = row == n_panels - 1

        ax_g.imshow(bar_rgb(gbar, gt_palette), aspect="auto", interpolation="nearest")
        ax_g.set_xticks([])
        ax_g.set_yticks([])
        ax_g.set_ylabel(name, fontsize=11, fontweight="bold", rotation=90,
                        labelpad=6, va="center")
        if row == 0:
            ax_g.set_title("GT", fontsize=8)

        if pbar is None:
            ax_p.axis("off")
        else:
            ax_p.imshow(bar_rgb(pbar, pred_palette), aspect="auto",
                        interpolation="nearest")
            ax_p.set_xticks([])
            ax_p.set_yticks([])
            n_clone = pd.Series(pbar).dropna().nunique()
            ax_p.set_title(f"n_clone={n_clone}", fontsize=8)

        if cnv is None:
            ax_h.text(0.5, 0.5, f"{name}: no output", ha="center", va="center",
                      fontsize=12, color="firebrick", transform=ax_h.transAxes)
            ax_h.set_xticks([])
            ax_h.set_yticks([])
        else:
            im = ax_h.imshow(np.ma.masked_invalid(cnv), aspect="auto",
                             interpolation="nearest", cmap=cmap,
                             vmin=-CLIP, vmax=CLIP)
            for s in sep_idx:
                ax_h.axvline(s - 0.5, color="black", lw=0.3)
            ax_h.set_yticks([])
            if is_bottom:
                ax_h.set_xticks([m for m, _ in chr_mids])
                ax_h.set_xticklabels(
                    [re.sub(r"^chr", "", c) for _, c in chr_mids], fontsize=7)
                ax_h.tick_params(axis="x", length=0)
                ax_h.set_xlabel("Genome (genes ordered by chromosome position)")
            else:
                ax_h.set_xticks([])

        for ax in (ax_g, ax_p, ax_h):
            if ax.axison:
                for r in grp_runs:
                    ax.axhline(r - 0.5, color="black", lw=0.5)

    if im is not None:
        cax = fig.add_axes([0.925, 0.35, 0.012, 0.3])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("CNV (log2 ratio, clipped ±1)", fontsize=9)

    def _sorted_handles(pal):
        return [Patch(facecolor=pal[c], edgecolor="black", lw=0.3, label=c)
                for c in sorted(pal,
                                key=lambda s: int(re.search(r"\d+", s).group()))]

    gt_handles = _sorted_handles(gt_palette)
    pred_handles = _sorted_handles(pred_palette)
    no_call = Patch(facecolor=str(NO_CALL_COLOR), edgecolor="black",
                    lw=0.3, label="no call")

    fig.legend(handles=gt_handles + [no_call], loc="lower center",
               ncol=min(12, len(gt_handles) + 1), fontsize=8, frameon=False,
               title="Ground truth subclone", bbox_to_anchor=(0.5, 0.015))
    fig.legend(handles=pred_handles + [no_call], loc="lower center",
               ncol=min(12, len(pred_handles) + 1), fontsize=8, frameon=False,
               title="Predicted subclone (per-tool clone indices)",
               bbox_to_anchor=(0.5, -0.04))

    fig.suptitle(
        f"Task 2 subclone -- ground truth vs tool CNV heatmaps  ({rep_name})\n"
        f"spots ordered by ground-truth subclone; left bars = GT | predicted; "
        f"InferCNV CNV = HMM17 TumorFreeNormal state",
        fontsize=12, fontweight="bold")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / "figS_task2C_cnv_heatmap"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {stem}.pdf / .png")


if __name__ == "__main__":
    main()
