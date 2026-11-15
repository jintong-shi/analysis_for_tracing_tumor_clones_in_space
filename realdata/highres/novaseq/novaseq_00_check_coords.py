#!/usr/bin/env python
"""Verify the spatial barcode -> (x, y) decoding for the Illumina NovaSeq X spatial data.

The download ships NO barcode coordinate file. The barcode IDs look like
    SBC:947855:4086433
and the hypothesis (inferred from range matching, not from vendor docs) is

    SBC:<y_nm>:<x_nm>     with the segmentation CSVs in micrometres  (nm = um * 1000)

This script tests that hypothesis per sample by rasterising the barcode point cloud
onto a coarse grid and overlaying the segmented-cell centroids. If the decoding is
right the two densities cover the same tissue outline.

Run:
  <CONDA_ROOT>/envs/bin2cell_env/bin/python \
    high_res_ST/novaseq/code/novaseq_00_check_coords.py
"""
import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from novaseq_samples import SAMPLES, OUTDIR

BIN_UM = 20.0


def barcode_xy(barcodes_path, max_n=None):
    """Decode SBC:<y_nm>:<x_nm> -> float arrays (x_um, y_um)."""
    ys, xs = [], []
    with gzip.open(barcodes_path, "rt") as fh:
        for i, line in enumerate(fh):
            if max_n is not None and i >= max_n:
                break
            _, a, b = line.rstrip("\n").split(":")
            ys.append(int(a))
            xs.append(int(b))
    return np.array(xs, dtype=np.float64) / 1000.0, np.array(ys, dtype=np.float64) / 1000.0


def hist2d(x, y, xlim, ylim, bin_um=BIN_UM):
    nx = int(np.ceil((xlim[1] - xlim[0]) / bin_um))
    ny = int(np.ceil((ylim[1] - ylim[0]) / bin_um))
    h, _, _ = np.histogram2d(
        x, y, bins=[nx, ny],
        range=[[xlim[0], xlim[0] + nx * bin_um], [ylim[0], ylim[0] + ny * bin_um]],
    )
    return h.T  # rows = y


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []
    fig, axes = plt.subplots(2, len(SAMPLES), figsize=(5 * len(SAMPLES), 11))

    for j, s in enumerate(SAMPLES):
        bx, by = barcode_xy(s["barcodes"])
        seg = pd.read_csv(s["cell_contours"])
        cent = seg.groupby("cell_id")[["vertex_x", "vertex_y"]].mean()

        xlim = (min(bx.min(), cent.vertex_x.min()), max(bx.max(), cent.vertex_x.max()))
        ylim = (min(by.min(), cent.vertex_y.min()), max(by.max(), cent.vertex_y.max()))

        hb = hist2d(bx, by, xlim, ylim)
        hc = hist2d(cent.vertex_x.values, cent.vertex_y.values, xlim, ylim)

        # overlap of occupied bins = the actual test
        ob, oc = hb > 0, hc > 0
        jacc = (ob & oc).sum() / (ob | oc).sum()
        cell_in_bc = (ob & oc).sum() / oc.sum()
        rows.append(dict(sample=s["name"], n_barcodes=len(bx), n_cells=len(cent),
                         bc_x=f"{bx.min():.0f}-{bx.max():.0f}", bc_y=f"{by.min():.0f}-{by.max():.0f}",
                         cell_x=f"{cent.vertex_x.min():.0f}-{cent.vertex_x.max():.0f}",
                         cell_y=f"{cent.vertex_y.min():.0f}-{cent.vertex_y.max():.0f}",
                         jaccard=round(jacc, 3), frac_cellbins_covered=round(cell_in_bc, 3)))

        ext = [xlim[0], xlim[0] + hb.shape[1] * BIN_UM, ylim[0] + hb.shape[0] * BIN_UM, ylim[0]]
        axes[0, j].imshow(np.log1p(hb), origin="upper", extent=ext, cmap="magma", aspect="equal")
        axes[0, j].set_title(f"{s['name']}\nbarcode density (log1p)", fontsize=9)
        axes[1, j].imshow(np.log1p(hc), origin="upper", extent=ext, cmap="viridis", aspect="equal")
        axes[1, j].set_title(f"segmented-cell centroid density\njaccard={jacc:.3f}", fontsize=9)
        for ax in (axes[0, j], axes[1, j]):
            ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
        print(rows[-1], flush=True)

    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "qc00_coord_decoding_check.png"), dpi=150)
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "qc00_coord_decoding_check.csv"), index=False)
    print("DONE ->", OUTDIR)


if __name__ == "__main__":
    main()
