#!/usr/bin/env python
"""Select normal / cancer epithelial cells in IDC_Grade3 for inferCNV.

A cell is kept when it is epithelial by clustering (major_celltype == "epi" at the
frozen res 0.3) AND falls inside a histology polygon:

    group = "normal"  -> inside one of the 4 "normal" polygons
    group = "cancer"  -> inside the "IDC" polygon

The geojson is in H&E pixels of <sample>_registered.ome.tiff; the OME metadata gives
PhysicalSizeX = PhysicalSizeY = 382 nm, so um = px * 0.382 puts the annotation in the
same frame as the cell centroids (x, y in um).

QC: nCount > 200 & nFeature > 150. Looser than the atera run (1000/800) because the
normal epithelium here is low-content (median 808 counts / 470 genes vs 3391 / 1697 in
the tumour) -- at the atera cutoff only 52 normal cells survive.

Env: bin2cell_env.
Run:
  MPLCONFIGDIR="$TMPDIR" <CONDA_ROOT>/envs/bin2cell_env/bin/python \
    high_res_ST/novaseq/code/novaseq_07_g3_select_regions.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path

SAMPLE = "IDC_Grade3"
PROJ = "<DATA_ROOT>/high_res_ST/novaseq"
GEOJSON = os.path.join(PROJ, "data", "annotation",
                       f"NovaSeqX_Human_Breast_{SAMPLE}_annotation.geojson")
CELLS = os.path.join(PROJ, "processed_data", "clustering", SAMPLE,
                     f"{SAMPLE}_annotated_cells_res0.3.csv")
OUT = os.path.join(PROJ, "processed_data", f"infercnv_{SAMPLE}")

PX_UM = 0.382          # PhysicalSizeX from the registered ome.tiff OME-XML
MIN_COUNT = 200
MIN_FEATURE = 150

os.makedirs(OUT, exist_ok=True)

## ---- 1. histology polygons -> um ---------------------------------
with open(GEOJSON) as fh:
    gj = json.load(fh)

cells = pd.read_csv(CELLS, index_col=0)
xy = cells[["x", "y"]].values
in_idc = np.zeros(len(cells), bool)
in_normal = np.zeros(len(cells), bool)

for feat in gj["features"]:
    rings = feat["geometry"]["coordinates"]
    assert feat["geometry"]["type"] == "Polygon" and len(rings) == 1, "unexpected geometry"
    inside = Path(np.asarray(rings[0]) * PX_UM).contains_points(xy)
    if feat["properties"]["name"] == "IDC":
        in_idc |= inside
    else:
        in_normal |= inside

cells["region"] = np.where(in_normal, "normal", np.where(in_idc, "IDC", "outside"))
print("cells per histology region:")
print(cells.region.value_counts())

## ---- 2. epithelial + QC ------------------------------------------
keep = (
    (cells.major_celltype == "epi")
    & (cells.nCount_spatial > MIN_COUNT)
    & (cells.nFeature_spatial > MIN_FEATURE)
    & cells.region.isin(["normal", "IDC"])
)
sel = cells.loc[keep].copy()
sel["group"] = np.where(sel.region == "normal", "normal", "cancer")

print(f"\nQC nCount > {MIN_COUNT} & nFeature > {MIN_FEATURE}, epithelial only:")
print(sel.group.value_counts())
print("\ngroup x cluster:")
print(pd.crosstab(sel.group, sel.cluster_frozen))

sel[["x", "y", "cluster_frozen", "nCount_spatial", "nFeature_spatial", "region", "group"]] \
    .to_csv(os.path.join(OUT, f"{SAMPLE}_epi_groups.csv"))

## ---- 3. selection map (visual check of the polygon mapping) ------
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(cells.x, cells.y, s=0.05, c="0.88", lw=0, rasterized=True)
ax.scatter(sel.x[sel.group == "cancer"], sel.y[sel.group == "cancer"],
           s=0.05, c="#B2182B", lw=0, rasterized=True, label=f"cancer epi (n={(sel.group == 'cancer').sum()})")
ax.scatter(sel.x[sel.group == "normal"], sel.y[sel.group == "normal"],
           s=6, c="#2166AC", lw=0, label=f"normal epi (n={(sel.group == 'normal').sum()})")
for feat in gj["features"]:
    ring = np.asarray(feat["geometry"]["coordinates"][0]) * PX_UM
    ax.plot(ring[:, 0], ring[:, 1], lw=0.8,
            c="k" if feat["properties"]["name"] == "IDC" else "#2166AC")
ax.invert_yaxis()
ax.set_aspect("equal")
ax.set_xlabel("x (um)")
ax.set_ylabel("y (um)")
ax.set_title(f"{SAMPLE}: epithelial cells by histology region")
ax.legend(markerscale=10, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, f"{SAMPLE}_epi_groups_map.png"), dpi=200)
print("\nwritten to", OUT)
