#!/usr/bin/env python3
"""Carrier-based ground-truth subclone extraction and QC plots.

Subclone rule:
  1. Start from tumor seed clone_id labels.
  2. For each evolution CNV event, restrict candidates to final tumor cells in
     the event cell's lineage.
  3. A candidate is a true carrier only when its final total-CN profile matches
     the event total copy number across enough affected genes.
  4. Only events whose true carriers reach --major-frac split subclones.

Copy-neutral allele-specific events are skipped: sc-genome_profile stores total
CN only.
"""

from __future__ import annotations

import argparse
import gc
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CN_VMAX = 5


def chrom_key(chrom: str) -> int:
    return CHROM_ORDER.index(chrom) if chrom in CHROM_ORDER else len(CHROM_ORDER)


def cell_num_to_name(cell_num: int) -> str:
    return f"Cell_{int(cell_num)}"


def cell_name_to_num(cell_name: str) -> int:
    return int(str(cell_name).replace("Cell_", ""))


def parse_cell_by_spot(path: Path) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for line in path.read_text().splitlines():
        if ":" not in line:
            continue
        spot, rhs = line.split(":", 1)
        vals = [v.strip() for v in rhs.split(",") if v.strip()]
        nums = []
        for val in vals:
            nums.append(cell_name_to_num(val))
        out[spot.strip()] = nums
    return out


def load_gene_positions(rep_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(rep_dir / "synthetic/genome.tsv", sep="\t", index_col=0)
    df = df[["gene", "chr", "start", "end", "mid"]].drop_duplicates("gene").copy()
    df["chr_order"] = df["chr"].map(chrom_key)
    return df.sort_values(["chr_order", "start"]).reset_index(drop=True)


def parse_event_genes(row: pd.Series, gene_pos: pd.DataFrame, profile_genes: set[str]) -> list[str]:
    genes: list[str] = []
    affected = row.get("affected_genes", "")
    if isinstance(affected, str) and affected and affected.lower() != "nan":
        genes.extend([g for g in affected.split(";") if g])

    genes = [g for g in genes if g in profile_genes]
    if genes:
        return list(dict.fromkeys(genes))

    chrom = str(row["Chromosome"])
    start = int(row["Start"])
    end = int(row["End"])
    mask = (
        gene_pos["chr"].eq(chrom)
        & (gene_pos["start"].astype(int) >= start)
        & (gene_pos["start"].astype(int) <= end)
    )
    pos_genes = gene_pos.loc[mask, "gene"].astype(str).tolist()
    return [g for g in pos_genes if g in profile_genes]


def load_sc_meta(rep_dir: Path) -> pd.DataFrame:
    sc = pd.read_csv(rep_dir / "synthetic/single_cell_data/sc-meta_data.tsv", sep="\t")
    sc = sc.rename(columns={sc.columns[0]: "cell_name"})
    sc["cell_num"] = sc["cell_name"].map(cell_name_to_num)
    sc["cell_type"] = sc["cell_type"].astype(str)
    sc["clone_id"] = pd.to_numeric(sc["clone_id"], errors="coerce").fillna(-1).astype(int)
    return sc


def load_lineage(rep_dir: Path) -> tuple[dict[int, int], dict[int, list[int]]]:
    path = rep_dir / "synthetic/lineage_all.tsv"
    if not path.exists():
        path = rep_dir / "synthetic/lineage.tsv"
    lin = pd.read_csv(path, sep="\t")
    # These tables are written with an unnamed row id, then child_id, parent_id.
    child_col = "child_id" if "child_id" in lin.columns else lin.columns[1]
    parent_col = "parent_id" if "parent_id" in lin.columns else lin.columns[2]
    lin = lin[[child_col, parent_col]].copy()
    lin.columns = ["child_id", "parent_id"]
    lin["child_id"] = pd.to_numeric(lin["child_id"], errors="coerce").astype("Int64")
    lin["parent_id"] = pd.to_numeric(lin["parent_id"], errors="coerce").astype("Int64")
    lin = lin.dropna(subset=["child_id"])

    parent_map = dict(zip(lin["child_id"].astype(int), lin["parent_id"].fillna(-1).astype(int)))
    children_map: dict[int, list[int]] = defaultdict(list)
    for child, parent in parent_map.items():
        if parent >= 0:
            children_map[parent].append(child)
    return parent_map, children_map


def descendants_of(founder: int, children_map: dict[int, list[int]], keep: set[int]) -> set[int]:
    out: set[int] = set()
    stack = [int(founder)]
    seen: set[int] = set()
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        if node in keep:
            out.add(node)
        stack.extend(children_map.get(node, []))
    return out


def event_label(row: pd.Series, used: set[str]) -> str:
    base = f"Evol_{int(row['cell_id'])}"
    label = base
    if label in used:
        label = f"{base}_row{int(row['row_id'])}"
    used.add(label)
    return label


def score_event_carriers(
    sc_profile_file: Path,
    genes: list[str],
    target_total_cn: int,
    candidate_cell_nums: set[int],
    profile_cell_names: set[str],
    carrier_gene_frac: float,
) -> tuple[set[int], dict[str, float]]:
    candidate_names = [
        cell_num_to_name(c) for c in sorted(candidate_cell_nums) if cell_num_to_name(c) in profile_cell_names
    ]
    if not genes or not candidate_names:
        return set(), {"mean_match_frac": 0.0, "max_match_frac": 0.0}

    df = pd.read_parquet(sc_profile_file, columns=["index"] + genes).set_index("index")
    sub = df.reindex(candidate_names)
    arr = sub.to_numpy(dtype=np.float32, copy=False)
    match = np.isclose(arr, float(target_total_cn), atol=0.01)
    frac = match.mean(axis=1)
    carriers = {
        cell_name_to_num(name)
        for name, f in zip(candidate_names, frac)
        if f >= carrier_gene_frac
    }
    stats = {
        "mean_match_frac": float(frac.mean()) if len(frac) else 0.0,
        "max_match_frac": float(frac.max()) if len(frac) else 0.0,
    }
    del df, sub, arr, match, frac
    gc.collect()
    return carriers, stats


def detect_event_carriers(
    rep_dir: Path,
    sc_meta: pd.DataFrame,
    gene_pos: pd.DataFrame,
    parent_map: dict[int, int],
    children_map: dict[int, list[int]],
    carrier_gene_frac: float,
) -> tuple[pd.DataFrame, dict[str, set[int]]]:
    event_log = pd.read_csv(rep_dir / "event_log.tsv", sep="\t")
    sc_profile_file = rep_dir / "synthetic/single_cell_data/sc-genome_profile.parquet"
    parquet_file = pq.ParquetFile(sc_profile_file)
    profile_cols = parquet_file.schema.names
    profile_genes = set(profile_cols) - {"index"}
    profile_cell_names = set(pd.read_parquet(sc_profile_file, columns=["index"])["index"].astype(str))

    tumor_cells = set(sc_meta.loc[sc_meta["cell_type"].str.lower().eq("tumor"), "cell_num"].astype(int))
    used_labels: set[str] = set()
    summary_rows: list[dict[str, object]] = []
    carriers_by_label: dict[str, set[int]] = {}

    for event_idx, row in event_log.iterrows():
        context = str(row.get("context", "")).lower()
        cell_id = int(row["cell_id"])
        target_total_cn = int(row["Major_Copy_Number"]) + int(row["Minor_Copy_Number"])
        label = event_label(row, used_labels) if context == "evolution" else f"{context}_{event_idx}"

        if context == "init_shared":
            candidates = set(tumor_cells)
        else:
            candidates = descendants_of(cell_id, children_map, tumor_cells)
            if not candidates and cell_id in tumor_cells:
                candidates = {cell_id}

        genes = parse_event_genes(row, gene_pos, profile_genes)
        status = "ok"
        carriers: set[int] = set()
        score_stats = {"mean_match_frac": 0.0, "max_match_frac": 0.0}

        if target_total_cn == 2:
            status = "skip_total_cn_neutral"
        elif not genes:
            status = "skip_no_profile_genes"
        elif not candidates:
            status = "skip_no_candidate_tumor_cells"
        else:
            carriers, score_stats = score_event_carriers(
                sc_profile_file=sc_profile_file,
                genes=genes,
                target_total_cn=target_total_cn,
                candidate_cell_nums=candidates,
                profile_cell_names=profile_cell_names,
                carrier_gene_frac=carrier_gene_frac,
            )
            if not carriers:
                status = "no_final_carriers"

        carriers_by_label[label] = carriers

        summary_rows.append(
            {
                "event_idx": event_idx,
                "label": label,
                "context": context,
                "iter": int(row["iter"]),
                "cell_id": cell_id,
                "row_id": int(row["row_id"]) if "row_id" in row and pd.notna(row["row_id"]) else np.nan,
                "chromosome": row["Chromosome"],
                "start": int(row["Start"]),
                "end": int(row["End"]),
                "target_total_cn": target_total_cn,
                "n_profile_genes": len(genes),
                "n_candidate_tumor_cells": len(candidates),
                "n_final_carrier_tumor_cells": len(carriers),
                "carrier_frac_all_tumor": len(carriers) / len(tumor_cells) if tumor_cells else 0.0,
                "mean_match_frac_candidate": score_stats["mean_match_frac"],
                "max_match_frac_candidate": score_stats["max_match_frac"],
                "status": status,
            }
        )

    return pd.DataFrame(summary_rows), carriers_by_label


def assign_subclones(
    sc_meta: pd.DataFrame,
    event_summary: pd.DataFrame,
    carriers_by_label: dict[str, set[int]],
    major_frac: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    tumor_cells = set(sc_meta.loc[sc_meta["cell_type"].str.lower().eq("tumor"), "cell_num"].astype(int))
    major = event_summary[
        (event_summary["context"].eq("evolution"))
        & (event_summary["status"].isin(["ok"]))
        & (event_summary["carrier_frac_all_tumor"] >= major_frac)
    ].copy()
    major_labels = major["label"].astype(str).tolist()

    rows: list[dict[str, object]] = []
    for rec in sc_meta.itertuples(index=False):
        cell_num = int(rec.cell_num)
        if str(rec.cell_type).lower() != "tumor":
            signature = "Normal"
            display_label = "Normal"
        else:
            carried = [label for label in major_labels if cell_num in carriers_by_label.get(label, set())]
            signature = f"Clone_{int(rec.clone_id)}|" + "|".join(carried)
            display_label = f"Clone_{int(rec.clone_id)}" if not carried else f"Clone_{int(rec.clone_id)}+" + "+".join(carried)
        rows.append(
            {
                "cell_name": rec.cell_name,
                "cell_num": cell_num,
                "x": float(rec.x),
                "y": float(rec.y),
                "cell_type": rec.cell_type,
                "clone_id": int(rec.clone_id),
                "subclone_signature": signature,
                "subclone_label": display_label,
            }
        )

    cell_df = pd.DataFrame(rows)
    tumor_signatures = sorted(cell_df.loc[cell_df["cell_type"].str.lower().eq("tumor"), "subclone_signature"].unique())
    sig_to_clone = {sig: f"clone{i + 1}" for i, sig in enumerate(tumor_signatures)}
    sig_to_clone["Normal"] = "Normal"
    cell_df["subclone"] = cell_df["subclone_signature"].map(sig_to_clone)
    cell_df["subclone_display"] = np.where(
        cell_df["subclone"].eq("Normal"),
        "Normal",
        cell_df["subclone"] + ": " + cell_df["subclone_label"],
    )

    major_out = major[
        [
            "label",
            "iter",
            "cell_id",
            "row_id",
            "chromosome",
            "start",
            "end",
            "target_total_cn",
            "n_profile_genes",
            "n_final_carrier_tumor_cells",
            "carrier_frac_all_tumor",
        ]
    ].copy()
    return cell_df, major_out, sig_to_clone


def build_spot_truth(rep_dir: Path, cell_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    spot_cells = parse_cell_by_spot(rep_dir / "synthetic/spatial_data/cell_by_spot.txt")
    cell_map = cell_df.set_index("cell_num")
    hard_rows: list[dict[str, object]] = []
    soft_rows: list[dict[str, object]] = []

    for spot_id, cell_nums in spot_cells.items():
        tumor_labels: list[str] = []
        tumor_display: list[str] = []
        for cell_num in cell_nums:
            if cell_num not in cell_map.index:
                continue
            rec = cell_map.loc[cell_num]
            if str(rec["cell_type"]).lower() == "tumor":
                tumor_labels.append(str(rec["subclone"]))
                tumor_display.append(str(rec["subclone_display"]))
        if not tumor_labels:
            continue

        counts = Counter(tumor_labels)
        total = sum(counts.values())
        major = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        display_by_clone = {}
        for clone, display in zip(tumor_labels, tumor_display):
            display_by_clone.setdefault(clone, display)
        hard_rows.append(
            {
                "spot_id": spot_id,
                "subclone": major,
                "subclone_display": display_by_clone.get(major, major),
            }
        )
        for clone, n in sorted(counts.items()):
            soft_rows.append(
                {
                    "spot_id": spot_id,
                    "subclone": clone,
                    "subclone_display": display_by_clone.get(clone, clone),
                    "n_cells": int(n),
                    "tumor_cells": int(total),
                    "prop": n / total if total else 0.0,
                }
            )

    return pd.DataFrame(hard_rows), pd.DataFrame(soft_rows)


def cnv_cmap_norm() -> tuple[ListedColormap, BoundaryNorm]:
    colors = ["#08306b", "#6baed6", "#f7f7f7", "#fb6a4a", "#cb181d", "#67000d"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, CN_VMAX + 1.5, 1), cmap.N)
    return cmap, norm


def add_chrom_separators(ax, gene_pos: pd.DataFrame) -> None:
    boundaries = []
    cur = gene_pos["chr"].iloc[0]
    boundaries.append((0, cur))
    for i, chrom in enumerate(gene_pos["chr"].values):
        if chrom != cur:
            boundaries.append((i, chrom))
            cur = chrom
    boundaries.append((len(gene_pos), None))
    for idx, _ in boundaries[1:-1]:
        ax.axvline(idx - 0.5, color="black", lw=0.4)
    label_pos = [((a + b) / 2, lab) for (a, lab), (b, _) in zip(boundaries[:-1], boundaries[1:])]
    ax.set_xticks([p for p, _ in label_pos])
    ax.set_xticklabels([re.sub(r"^chr", "", lab) for _, lab in label_pos], fontsize=7)
    ax.tick_params(axis="x", which="both", length=0, pad=2)


def runs_of(labels: pd.Series) -> list[tuple[int, int, str]]:
    labels = labels.reset_index(drop=True)
    runs = []
    start = 0
    for i in range(1, len(labels)):
        if labels.iat[i] != labels.iat[i - 1]:
            runs.append((start, i, labels.iat[i - 1]))
            start = i
    runs.append((start, len(labels), labels.iat[-1]))
    return runs


def stable_colors(labels: list[str]) -> dict[str, tuple]:
    palette = plt.get_cmap("tab10")(np.linspace(0, 1, max(10, len(labels))))
    return {label: tuple(palette[i % len(palette)]) for i, label in enumerate(labels)}


def plot_cnv_heatmap(
    matrix: np.ndarray,
    gene_pos: pd.DataFrame,
    row_groups: pd.Series,
    title: str,
    ylabel: str,
    out_path: Path,
) -> None:
    cmap, norm = cnv_cmap_norm()
    matrix = np.clip(matrix, 0, CN_VMAX)
    n_rows = matrix.shape[0]
    fig = plt.figure(figsize=(16, max(5.0, min(11, n_rows / 220))))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[0.18, 14, 0.25], wspace=0.04)
    ax_side = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1], sharey=ax_side)
    cax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    add_chrom_separators(ax, gene_pos)
    row_runs = runs_of(row_groups)
    for _, end, _ in row_runs[:-1]:
        ax.axhline(end - 0.5, color="black", lw=0.4)
    ax.set_yticks([])
    ax.set_xlabel("Gene (genomic position)")
    ax.set_title(title)

    labels = list(dict.fromkeys([r[2] for r in row_runs]))
    color_map = stable_colors(labels)
    side = np.zeros((n_rows, 1, 3), dtype=float)
    for start, end, label in row_runs:
        side[start:end, 0, :] = color_map[label][:3]
    ax_side.imshow(side, aspect="auto", interpolation="nearest")
    for _, end, _ in row_runs[:-1]:
        ax_side.axhline(end - 0.5, color="black", lw=0.4)
    ax_side.set_xticks([])
    ax_side.set_yticks([])
    ax_side.set_ylabel(ylabel)

    cbar = fig.colorbar(im, cax=cax, ticks=range(0, CN_VMAX + 1))
    cbar.set_ticklabels([str(i) for i in range(CN_VMAX)] + [f"{CN_VMAX}+"])
    cbar.set_label("Copy number")
    handles = [Patch(facecolor=color_map[label], edgecolor="black", lw=0.3, label=label) for label in labels]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=min(5, len(handles)),
        fontsize=8,
        frameon=False,
        title="Subclone",
    )
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sc_heatmap(rep_dir: Path, out_dir: Path, cell_df: pd.DataFrame, gene_pos: pd.DataFrame, max_cells: int) -> None:
    tumor = cell_df[cell_df["cell_type"].str.lower().eq("tumor")].copy()
    if max_cells and len(tumor) > max_cells:
        rng = np.random.default_rng(0)
        keep_idx = []
        for _, sub in tumor.groupby("subclone", sort=False):
            n = max(1, int(round(len(sub) / len(tumor) * max_cells)))
            n = min(n, len(sub))
            keep_idx.extend(rng.choice(sub.index.values, size=n, replace=False))
        tumor = tumor.loc[keep_idx].copy()
    tumor = tumor.sort_values(["subclone", "cell_num"], kind="stable").reset_index(drop=True)

    gene_order = gene_pos["gene"].astype(str).tolist()
    sc_geno = pd.read_parquet(rep_dir / "synthetic/single_cell_data/sc-genome_profile.parquet")
    sc_geno = sc_geno.set_index("index").reindex(columns=gene_order)
    mat = sc_geno.loc[tumor["cell_name"].tolist()].to_numpy()
    plot_cnv_heatmap(
        mat,
        gene_pos,
        tumor["subclone"],
        f"Ground-truth CNV (single-cell) - {rep_dir.name}",
        f"Cells (n={len(tumor)})",
        out_dir / "groundtruth_sc_heatmap.png",
    )
    del sc_geno, mat
    gc.collect()


def plot_spot_heatmap(rep_dir: Path, out_dir: Path, hard_df: pd.DataFrame, gene_pos: pd.DataFrame) -> None:
    gene_order = gene_pos["gene"].astype(str).tolist()
    st_geno = pd.read_parquet(rep_dir / "synthetic/spatial_data/st-genome_profile.parquet")
    st_geno = st_geno.set_index("index").reindex(index=gene_order).T
    hard = hard_df.set_index("spot_id")
    spots = [s for s in st_geno.index if s in hard.index]
    order = pd.DataFrame({"spot_id": spots, "subclone": hard.loc[spots, "subclone"].values})
    order = order.sort_values(["subclone", "spot_id"], kind="stable").reset_index(drop=True)
    mat = np.round(st_geno.loc[order["spot_id"].tolist()].to_numpy()).astype(int)
    plot_cnv_heatmap(
        mat,
        gene_pos,
        order["subclone"],
        f"Ground-truth CNV (spot, rounded) - {rep_dir.name}",
        f"Spots (n={len(order)})",
        out_dir / "groundtruth_spot_heatmap.png",
    )
    del st_geno, mat
    gc.collect()


def plot_subclone_spatial(rep_dir: Path, out_dir: Path, cell_df: pd.DataFrame, n_subclones: int) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    normal_mask = ~cell_df["cell_type"].str.lower().eq("tumor").to_numpy()
    ax.scatter(
        cell_df.loc[normal_mask, "x"],
        cell_df.loc[normal_mask, "y"],
        c="#e0e0e0",
        s=6,
        alpha=0.3,
        rasterized=True,
        linewidths=0,
    )
    tumor = cell_df[~normal_mask].copy()
    labels = sorted(tumor["subclone"].unique())
    color_map = stable_colors(labels)
    for label in labels:
        sub = tumor[tumor["subclone"].eq(label)]
        ax.scatter(sub["x"], sub["y"], c=[color_map[label]], s=10, alpha=0.6, rasterized=True, linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"{rep_dir.name}: Subclone Spatial ({n_subclones} subclones)", fontsize=13, fontweight="bold")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#e0e0e0", markersize=8, label="Normal")]
    handles.extend(
        Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[label], markersize=8, label=label)
        for label in labels
    )
    ax.legend(handles=handles, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.tight_layout()
    fig.savefig(out_dir / "subclone_spatial.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_cnv_events_spatial(rep_dir: Path, out_dir: Path, cell_df: pd.DataFrame, event_summary: pd.DataFrame, carriers_by_label: dict[str, set[int]]) -> None:
    counts = Counter()
    for label, carriers in carriers_by_label.items():
        status = event_summary.loc[event_summary["label"].eq(label), "status"]
        if not status.empty and status.iloc[0] == "ok":
            for cell_num in carriers:
                counts[cell_num] += 1
    values = cell_df["cell_num"].map(lambda x: counts.get(int(x), 0)).to_numpy()

    fig, ax = plt.subplots(figsize=(8, 7))
    normal_mask = ~cell_df["cell_type"].str.lower().eq("tumor").to_numpy()
    ax.scatter(cell_df.loc[normal_mask, "x"], cell_df.loc[normal_mask, "y"], c="#e0e0e0", s=6, alpha=0.3, rasterized=True, linewidths=0)
    tumor_mask = ~normal_mask
    if tumor_mask.any():
        sc = ax.scatter(
            cell_df.loc[tumor_mask, "x"],
            cell_df.loc[tumor_mask, "y"],
            c=values[tumor_mask],
            cmap="YlOrRd",
            s=10,
            alpha=0.6,
            rasterized=True,
            linewidths=0,
            vmin=0,
        )
        plt.colorbar(sc, ax=ax, shrink=0.7, label="# visible total-CN CNV events")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"{rep_dir.name}: Total Visible CNV Events", fontsize=13, fontweight="bold")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.tight_layout()
    fig.savefig(out_dir / "cnv_events_spatial.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_spot_major_subclone(rep_dir: Path, out_dir: Path, hard_df: pd.DataFrame) -> None:
    st_meta = pd.read_parquet(rep_dir / "synthetic/spatial_data/st-meta.parquet").rename(columns={"index": "spot_id"})
    labels = st_meta[["spot_id", "x", "y"]].merge(hard_df[["spot_id", "subclone"]], on="spot_id", how="left")
    labels["subclone"] = labels["subclone"].fillna("Normal")
    unique = sorted(labels["subclone"].unique())
    tumor_labels = [x for x in unique if x not in {"Normal", "Empty", "Unknown"}]
    color_map = stable_colors(tumor_labels)
    color_map.update({"Normal": (0.85, 0.85, 0.85, 1.0), "Empty": (1.0, 1.0, 1.0, 1.0), "Unknown": (0.95, 0.95, 0.95, 1.0)})

    fig, ax = plt.subplots(figsize=(8, 7))
    for label in unique:
        sub = labels[labels["subclone"].eq(label)]
        ax.scatter(sub["x"], sub["y"], c=[color_map.get(label, (0.5, 0.5, 0.5, 1.0))], s=120, edgecolors="grey", linewidths=0.5)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"{rep_dir.name}: Major Subclone per Spot", fontsize=13, fontweight="bold")
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map.get(label, (0.5, 0.5, 0.5, 1.0)), markeredgecolor="grey", markersize=10, label=label)
        for label in unique
    ]
    ax.legend(handles=handles, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.tight_layout()
    fig.savefig(out_dir / "st_major_subclone_spatial.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_spot_subclone_overlay(rep_dir: Path, out_dir: Path, soft_df: pd.DataFrame) -> None:
    st_meta = pd.read_parquet(rep_dir / "synthetic/spatial_data/st-meta.parquet").rename(columns={"index": "spot_id"})
    labels = sorted(soft_df["subclone"].unique()) if len(soft_df) else []
    color_map = stable_colors(labels)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(st_meta["x"], st_meta["y"], c="#eeeeee", s=115, edgecolors="grey", linewidths=0.4)
    for label in labels:
        cur = soft_df[soft_df["subclone"].eq(label)].merge(st_meta[["spot_id", "x", "y"]], on="spot_id", how="left")
        ax.scatter(cur["x"], cur["y"], c=[color_map[label]], s=30 + 90 * cur["prop"], alpha=0.65, linewidths=0, label=label)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"{rep_dir.name}: Spot Subclone Composition", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.tight_layout()
    fig.savefig(out_dir / "st_subclone_count_spatial.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_n_clones(out_dir: Path, rep_name: str, n_clones: int, major_frac: float, event_summary: pd.DataFrame) -> None:
    n_major = int(((event_summary["context"].eq("evolution")) & (event_summary["status"].eq("ok")) & (event_summary["carrier_frac_all_tumor"] >= major_frac)).sum())
    lines = [
        "# STARCH n_clones from carrier-based ground truth",
        f"rep={rep_name}",
        f"n_clones={n_clones}",
        f"major_subclone_min_cell_frac={major_frac:.2f}",
        f"n_major_evolution_events={n_major}",
        "status=ok",
    ]
    (out_dir / "n_clones.txt").write_text("\n".join(lines) + "\n")


def write_production_truth(rep_dir: Path, hard_df: pd.DataFrame, soft_df: pd.DataFrame) -> tuple[Path, Path]:
    spatial_dir = rep_dir / "synthetic" / "spatial_data"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    hard_out = spatial_dir / "subclone_groundtruth_spot_hard.csv"
    soft_out = spatial_dir / "subclone_groundtruth_spot_soft_long.csv"
    hard_df.to_csv(hard_out, index=False)
    soft_df.to_csv(soft_out, index=False)
    return hard_out, soft_out


def run(
    rep_dir: Path,
    out_root: Path,
    major_frac: float,
    carrier_gene_frac: float,
    max_cells: int,
    write_production: bool = False,
) -> dict[str, object]:
    rep_dir = rep_dir.resolve()
    out_dir = out_root / rep_dir.name
    qc_dir = out_dir / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    gene_pos = load_gene_positions(rep_dir)
    sc_meta = load_sc_meta(rep_dir)
    parent_map, children_map = load_lineage(rep_dir)
    event_summary, carriers_by_label = detect_event_carriers(
        rep_dir=rep_dir,
        sc_meta=sc_meta,
        gene_pos=gene_pos,
        parent_map=parent_map,
        children_map=children_map,
        carrier_gene_frac=carrier_gene_frac,
    )
    cell_df, major_events, sig_to_clone = assign_subclones(sc_meta, event_summary, carriers_by_label, major_frac)
    hard_df, soft_df = build_spot_truth(rep_dir, cell_df)

    n_clones = int(cell_df.loc[cell_df["cell_type"].str.lower().eq("tumor"), "subclone"].nunique())
    event_summary.to_csv(out_dir / "event_carrier_summary.csv", index=False)
    major_events.to_csv(out_dir / "major_events.csv", index=False)
    pd.DataFrame({"subclone_signature": list(sig_to_clone.keys()), "subclone": list(sig_to_clone.values())}).to_csv(
        out_dir / "subclone_signature_map.csv", index=False
    )
    cell_df.to_csv(out_dir / "cell_subclone.csv", index=False)
    hard_df.to_csv(out_dir / "subclone_groundtruth_spot_hard.csv", index=False)
    soft_df.to_csv(out_dir / "subclone_groundtruth_spot_soft_long.csv", index=False)
    write_n_clones(qc_dir, rep_dir.name, n_clones, major_frac, event_summary)
    hard_prod_out = ""
    soft_prod_out = ""
    if write_production:
        hard_prod_out, soft_prod_out = write_production_truth(rep_dir, hard_df, soft_df)

    plot_sc_heatmap(rep_dir, qc_dir, cell_df, gene_pos, max_cells)
    plot_spot_heatmap(rep_dir, qc_dir, hard_df, gene_pos)
    plot_cnv_events_spatial(rep_dir, qc_dir, cell_df, event_summary, carriers_by_label)
    plot_subclone_spatial(rep_dir, qc_dir, cell_df, n_clones)
    plot_spot_major_subclone(rep_dir, qc_dir, hard_df)
    plot_spot_subclone_overlay(rep_dir, qc_dir, soft_df)

    print(f"[gt] rep={rep_dir.name}")
    print(f"[gt] n_clones={n_clones}")
    print(f"[gt] output={out_dir}")
    print("[gt] major events:")
    if len(major_events):
        print(major_events.to_string(index=False))
    else:
        print("  none")
    return {
        "rep": rep_dir.name,
        "out_dir": str(out_dir),
        "hard_out": str(hard_prod_out),
        "soft_out": str(soft_prod_out),
        "n_spots_with_subclone": int(len(hard_df)),
        "n_subclones_cell_level": n_clones,
        "n_subclones": int(hard_df["subclone"].nunique()) if len(hard_df) else 0,
        "n_tumor_cells": int(cell_df["cell_type"].str.lower().eq("tumor").sum()),
        "n_major_evolution_events": int(len(major_events)),
        "major_subclone_min_cell_frac": float(major_frac),
        "status": "ok",
    }


def run_root(
    root_dir: Path,
    out_root: Path,
    major_frac: float,
    carrier_gene_frac: float,
    max_cells: int,
    write_production: bool,
) -> None:
    rep_dirs = sorted(root_dir.glob("Rep_*"), key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)])
    if not rep_dirs:
        raise SystemExit(f"No Rep_* folders found under: {root_dir}")

    rows = []
    for rep_dir in rep_dirs:
        rows.append(
            run(
                rep_dir=rep_dir,
                out_root=out_root,
                major_frac=major_frac,
                carrier_gene_frac=carrier_gene_frac,
                max_cells=max_cells,
                write_production=write_production,
            )
        )

    results_dir = root_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(results_dir / "subclone_groundtruth_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rep-dir", type=Path)
    parser.add_argument("--root-dir", type=Path)
    parser.add_argument("--out-root", type=Path)
    parser.add_argument("--major-frac", type=float, default=0.05)
    parser.add_argument("--carrier-gene-frac", type=float, default=0.80)
    parser.add_argument("--max-cells", type=int, default=3000)
    parser.add_argument(
        "--write-production",
        action="store_true",
        help="Also write subclone_groundtruth_spot_*.csv into each rep's synthetic/spatial_data.",
    )
    args = parser.parse_args()
    if (args.rep_dir is None) == (args.root_dir is None):
        raise SystemExit("Provide exactly one of --rep-dir or --root-dir.")

    if args.rep_dir is not None:
        if args.out_root is None:
            raise SystemExit("--out-root is required with --rep-dir.")
        run(
            args.rep_dir,
            args.out_root.resolve(),
            args.major_frac,
            args.carrier_gene_frac,
            args.max_cells,
            write_production=args.write_production,
        )
    else:
        root_dir = args.root_dir.resolve()
        out_root = args.out_root.resolve() if args.out_root is not None else root_dir
        run_root(
            root_dir=root_dir,
            out_root=out_root,
            major_frac=args.major_frac,
            carrier_gene_frac=args.carrier_gene_frac,
            max_cells=args.max_cells,
            write_production=True if args.out_root is None else args.write_production,
        )


if __name__ == "__main__":
    main()
