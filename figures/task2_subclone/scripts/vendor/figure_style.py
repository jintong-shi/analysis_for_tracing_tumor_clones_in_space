"""Shared style constants and helpers for spatial-CNV benchmark figures.

Vendored copy of the pre-release `scripts/figures/figure_style.py`, kept only
because `figS_task2C_prep.py` imports the original figS_task2C loaders verbatim
(see `vendor/subclone_heatmap.py`). Not a public API: figure scripts in this
repo inline the style helpers they need. Only the module paths were changed --
`figure_settings.tsv` now resolves to the repo's canonical `figures/utils/`
copy (a superset of the pre-release table; no row this module reads differs),
and `project_path()` is rooted at `<DATA_ROOT>` (override at run time with the
`DATA_ROOT` environment variable instead of substituting the placeholder).
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FIGURE_SHARED_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get("DATA_ROOT", "<DATA_ROOT>"))
FIGURE_SETTINGS_PATH = FIGURE_SHARED_DIR.parents[2] / "utils" / "figure_settings.tsv"


def project_path(*parts: str) -> Path:
    return PROJECT_ROOT.joinpath(*parts)


def _decode_label(value: str | None) -> str | None:
    if value is None:
        return None
    return value.replace("\\n", "\n")


def load_figure_settings() -> list[dict[str, str | None]]:
    with FIGURE_SETTINGS_PATH.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key in ("label", "alias"):
            row[key] = _decode_label(row.get(key))
        if row.get("color") == "":
            row["color"] = None
    return rows


FIGURE_SETTINGS = load_figure_settings()


def figure_settings(setting_type: str | None = None,
                    scope: str | None = None) -> list[dict[str, str | None]]:
    rows = FIGURE_SETTINGS
    if setting_type is not None:
        rows = [row for row in rows if row["setting_type"] == setting_type]
    if scope is not None:
        rows = [row for row in rows if row["scope"] == scope]
    return sorted(rows, key=lambda row: (float(row["order"] or 1e9), row["id"]))


def setting_order(scope: str, setting_type: str = "tool_order",
                  value: str = "id") -> list[str]:
    return [str(row[value]) for row in figure_settings(setting_type, scope)]


TOOL_TABLE = figure_settings("tool", "global")
TOOL_DISPLAY = {str(row["id"]): str(row["label"]) for row in TOOL_TABLE}
TOOL_COLORS = {str(row["id"]): str(row["color"]) for row in TOOL_TABLE}


def _tool_aliases() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for row in TOOL_TABLE:
        tool_id = str(row["id"])
        values = [tool_id, str(row["label"])]
        alias = row.get("alias")
        if alias:
            values.extend(str(alias).split("|"))
        for value in values:
            value = value.strip()
            if value:
                aliases.setdefault(value, tool_id)
    return aliases


TOOL_ALIAS_TO_ID = _tool_aliases()
TOOL_DISPLAY_BY_ALIAS = {
    alias: TOOL_DISPLAY[tool_id] for alias, tool_id in TOOL_ALIAS_TO_ID.items()
}
TOOL_COLORS_BY_ALIAS = {
    alias: TOOL_COLORS[tool_id] for alias, tool_id in TOOL_ALIAS_TO_ID.items()
}


def tool_id(value: str) -> str:
    return TOOL_ALIAS_TO_ID.get(str(value), str(value))


def tool_label(value: str) -> str:
    return TOOL_DISPLAY.get(tool_id(value), str(value))


def tool_color(value: str) -> str | None:
    return TOOL_COLORS.get(tool_id(value))


TOOL_ORDER_REAL_TASK3 = setting_order("real_task3")
TOOL_ORDER_SYNTH_TASK2_ID = setting_order("synth_task2_summary")
TOOL_ORDER_SYNTH_TASK2 = [tool_label(value) for value in TOOL_ORDER_SYNTH_TASK2_ID]
TOOL_ORDER_SYNTH_TASK2_OVERVIEW_ID = setting_order("synth_task2_overview_ggplot")
TOOL_ORDER_SYNTH_TASK2_OVERVIEW = [
    tool_label(value) for value in TOOL_ORDER_SYNTH_TASK2_OVERVIEW_ID
]
TOOL_ORDER_SYNTH_TASK2_HEATMAP_PY_ID = setting_order("synth_task2_heatmap_py")
TOOL_ORDER_SYNTH_TASK2_HEATMAP_PY = [
    tool_label(value) for value in TOOL_ORDER_SYNTH_TASK2_HEATMAP_PY_ID
]
TOOL_COLORS_SYNTH_TASK2 = {
    label: str(tool_color(label)) for label in TOOL_ORDER_SYNTH_TASK2
}

# Backwards-compatible name for the current Python heatmap script.
TOOLS = TOOL_ORDER_SYNTH_TASK2_HEATMAP_PY

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CNV_CLIP = 1.0
CLIP = CNV_CLIP

SUBCLONE_GT_CMAP = "Set1"
SUBCLONE_PRED_CMAP = "tab20"
CNV_CMAP = "RdBu_r"
NO_CALL_COLOR = 0.85

DEFAULT_REP = str(project_path(
    "synthetic",
    "processed_data",
    "random_batch_BAF",
    "generate_code",
    "Complete_Batch_v2",
    "Rep_13_Clones_3_Shared_3",
))
DEFAULT_OUT = str(project_path(
    "synthetic",
    "results",
    "task2_subclone",
    "figures",
))


def chrom_key(chrom: str) -> int:
    return CHROM_ORDER.index(chrom) if chrom in CHROM_ORDER else len(CHROM_ORDER)


def _sort_clones(labels: Iterable[object]) -> list[str]:
    return sorted(
        {str(value) for value in labels if pd.notna(value)},
        key=lambda value: (
            int(match.group(1)) if (match := re.search(r"(\d+)", value)) else 1e9,
            value,
        ),
    )


def clone_palette(labels: Iterable[object], cmap_name: str) -> dict[str, object]:
    """Map clone labels to colors from a qualitative matplotlib colormap."""
    clones = _sort_clones(labels)
    cmap = plt.get_cmap(cmap_name)
    n_colors = getattr(cmap, "N", 20)
    return {clone: cmap(index % n_colors) for index, clone in enumerate(clones)}


def bar_rgb(labels: Iterable[object], palette: dict[str, object],
            missing: float = NO_CALL_COLOR) -> np.ndarray:
    """Return an (n, 1, 3) RGB strip; missing labels are light grey."""
    labels = list(labels)
    rgb = np.full((len(labels), 1, 3), missing, dtype=float)
    for index, label in enumerate(labels):
        label_text = str(label)
        if pd.notna(label) and label_text in palette:
            rgb[index, 0, :] = palette[label_text][:3]
    return rgb


def chrom_boundaries(gene_pos: pd.DataFrame):
    """Return chromosome boundary indices and midpoint labels."""
    chroms = gene_pos["chr"].to_numpy()
    bounds = [0] + [
        index for index in range(1, len(chroms))
        if chroms[index] != chroms[index - 1]
    ]
    bounds.append(len(chroms))
    mids = [
        ((bounds[index] + bounds[index + 1]) / 2, chroms[bounds[index]])
        for index in range(len(bounds) - 1)
    ]
    return bounds[1:-1], mids
