#!/usr/bin/env python
"""Task 2 subclone -- ground-truth vs tool CNV heatmaps.

Vendored, unmodified except for paths: this is the pre-release
`paper/figures/fig3_syn_task2_subclone_ari/subclone_heatmap.py`, kept here only
because `figS_task2C_prep.py` imports its raw loaders and panel table verbatim
so the compact figure stays bit-identical to the frozen original. Do not use it
as the figure script -- `figS_task2C_cnv_heatmap.py` is the released plotter.
Changed here: `figure_style` is imported from the sibling vendored copy, and the
data roots resolve under `<DATA_ROOT>/synthetic/`.

Combined figure: a ground-truth CNV heatmap on top, then each tool's predicted
CNV heatmap stacked below. Every panel shares the same spot ordering (spots
ordered by ground-truth subclone) and the same genomic gene order on x, so the
ground truth can be compared against each tool at a glance.

Adapted for the numbat_dirichlet_fix_test tools_batch layout:
  * ground truth (genome.tsv, st-genome_profile.parquet) lives in the SOURCE rep
    runs/<rep>/synthetic, not under tools_batch -> two roots (--rep_dir tool
    outputs, --syn_dir synthetic ground truth).
  * GT spot labels come from the co-located subclone_groundtruth_spot_hard.csv.
  * InferCNV CNV panel uses the HMM17_TumorFreeNormal *state* file
    (processed_state_spots_x_genes.csv; loss/neutral/gain -> -1/0/+1 on the same
    RdBu scale). Its predicted-subclone bar still comes from the noHMM run.

Style constants (tool order, colormaps) come from scripts/figures/figure_style.py.

Usage:
  python subclone_heatmap.py [--rep_dir DIR] [--syn_dir DIR] [--out_dir DIR]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
# Arial per FIGURE_STANDARDS; Nimbus Sans is the metric-compatible fallback on
# machines without Arial (e.g. st-a100). matplotlib falls through the list.
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Nimbus Sans", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import (
    CLIP,
    CNV_CMAP,
    NO_CALL_COLOR,
    SUBCLONE_GT_CMAP,
    SUBCLONE_PRED_CMAP,
    bar_rgb,
    chrom_boundaries,
    chrom_key,
    clone_palette,
    project_path,
)

# Per-panel table (one row per tool below the ground-truth panel). Columns:
#   display name | CNV kind ("log2" | "state") | CNV subdir | subclone subdir.
# The InferCNV row takes its predicted-subclone bar from the noHMM run (its only
# subclone labels) and its CNV from the HMM17 run.
PANELS = [
    ("InferCNV (HMM state)", "state", "InferCNV_HMM17_TumorFreeNormal", "InferCNV_noHMM_TumorFreeNormal"),
    ("SCEVAN",                     "log2",  "SCEVAN_TumorFreeNormal",         "SCEVAN_TumorFreeNormal"),
    ("STARCH",                     "log2",  "STARCH",                         "STARCH"),
    ("Numbat",                     "log2",  "Numbat_TumorFreeNormal",         "Numbat_TumorFreeNormal"),
    ("CalicoST",                   "log2",  "CalicoST",                       "CalicoST"),
]
STATE_MAP = {"loss": -CLIP, "neutral": 0.0, "gain": CLIP}

DEFAULT_REP = str(project_path(
    "synthetic", "numbat_dirichlet_fix_test", "runs",
    "tools_batch", "Rep_13_Clones_3_Shared_3", "Rep_13_Clones_3_Shared_3"))
DEFAULT_SYN = str(project_path(
    "synthetic", "numbat_dirichlet_fix_test", "runs",
    "Rep_13_Clones_3_Shared_3", "synthetic"))
DEFAULT_OUT = str(Path(__file__).resolve().parent)


def load_gene_order(syn_dir: Path) -> pd.DataFrame:
    """Genes ordered by genomic position; columns [gene, chr, start]."""
    g = pd.read_csv(syn_dir / "genome.tsv", sep="\t")
    g = g[["gene", "chr", "start"]].drop_duplicates(subset="gene", keep="first")
    g["chr_order"] = g["chr"].map(chrom_key)
    return g.sort_values(["chr_order", "start"]).reset_index(drop=True)


def load_gt_cnv(syn_dir: Path, genes: list[str], spots: list[str]) -> np.ndarray:
    """Ground-truth CNV as log2(CN/2): spots x genes, clipped to [-CLIP, CLIP]."""
    st = pd.read_parquet(syn_dir / "spatial_data" / "st-genome_profile.parquet")
    st = st.set_index(st.columns[0])
    mat = st.reindex(index=genes, columns=spots).to_numpy(dtype=float).T
    with np.errstate(divide="ignore", invalid="ignore"):
        log2 = np.log2(np.clip(mat, 1e-3, None) / 2.0)
    return np.clip(log2, -CLIP, CLIP)


def load_tool_cnv(rep_dir: Path, subdir: str, genes: list[str],
                  spots: list[str]) -> np.ndarray | None:
    """Tool predicted CNV (already log2): spots x genes, clipped. None if absent."""
    f = rep_dir / "Tool_Outputs" / subdir / "processed_log2_spots_x_genes.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f, index_col=0)
    mat = d.reindex(index=spots, columns=genes).to_numpy(dtype=float)
    # STARCH writes a large negative sentinel (~-39.86) for undefined entries;
    # treat anything below a real log2 range as missing rather than a deletion.
    mat[mat < -10.0] = np.nan
    return np.clip(mat, -CLIP, CLIP)


def load_infercnv_state(rep_dir: Path, subdir: str, genes: list[str],
                        spots: list[str]) -> np.ndarray | None:
    """InferCNV HMM state (loss/neutral/gain) -> {-CLIP, 0, +CLIP}: spots x genes."""
    f = rep_dir / "Tool_Outputs" / subdir / "processed_state_spots_x_genes.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f, index_col=0).reindex(index=spots, columns=genes)
    arr = d.to_numpy(dtype=object)
    out = np.full(arr.shape, np.nan, dtype=float)
    for key, val in STATE_MAP.items():
        out[arr == key] = val
    return out


def load_pred_labels(rep_dir: Path, subdir: str, spots: list[str]) -> pd.Series:
    """Predicted subclone label per spot (NaN where the tool gave no call)."""
    f = rep_dir / "Tool_Outputs" / subdir / "subclone_spots_x_subclones.csv"
    if not f.exists():
        return pd.Series(index=spots, dtype=object)
    p = pd.read_csv(f).set_index("spot_id")["subclone"]
    return p.reindex(spots)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep_dir", default=DEFAULT_REP, help="tool-outputs root")
    ap.add_argument("--syn_dir", default=DEFAULT_SYN, help="synthetic GT root")
    ap.add_argument("--out_dir", default=DEFAULT_OUT)
    args = ap.parse_args()

    rep_dir = Path(args.rep_dir)
    syn_dir = Path(args.syn_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_name = rep_dir.name

    gene_pos = load_gene_order(syn_dir)
    genes = gene_pos["gene"].tolist()

    gt_hard = pd.read_csv(rep_dir / "subclone_groundtruth_spot_hard.csv")
    gt_hard = gt_hard.sort_values(["subclone", "spot_id"], kind="stable")
    spots = gt_hard["spot_id"].tolist()
    gt_lab = gt_hard["subclone"].tolist()

    panels = [("Ground truth", load_gt_cnv(syn_dir, genes, spots), gt_lab, None)]
    pred_all = []
    for name, kind, cnv_dir, sub_dir in PANELS:
        if kind == "state":
            cnv = load_infercnv_state(rep_dir, cnv_dir, genes, spots)
        else:
            cnv = load_tool_cnv(rep_dir, cnv_dir, genes, spots)
        pred = load_pred_labels(rep_dir, sub_dir, spots)
        pred_all.append(pred)
        panels.append((name, cnv, gt_lab, pred))

    gt_palette = clone_palette(gt_lab, SUBCLONE_GT_CMAP)
    pred_palette = clone_palette(
        [x for p in pred_all for x in p], SUBCLONE_PRED_CMAP
    )

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
               title="Ground truth subclone",
               bbox_to_anchor=(0.5, 0.015))
    fig.legend(handles=pred_handles + [no_call], loc="lower center",
               ncol=min(12, len(pred_handles) + 1), fontsize=8, frameon=False,
               title="Predicted subclone (per-tool clone indices)",
               bbox_to_anchor=(0.5, -0.04))

    fig.suptitle(
        f"Task 2 subclone -- ground truth vs tool CNV heatmaps  ({rep_name})\n"
        f"spots ordered by ground-truth subclone; left bars = GT | predicted; "
        f"InferCNV CNV = HMM17 TumorFreeNormal state",
        fontsize=12, fontweight="bold")

    stem = out_dir / f"subclone_cnv_heatmap_{rep_name}"
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {stem}.pdf / .png")


if __name__ == "__main__":
    main()
