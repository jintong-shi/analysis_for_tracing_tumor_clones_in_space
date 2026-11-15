#!/usr/bin/env python
"""Deliverable 1a — build a cell x gene matrix for each NovaSeq X Illumina spatial sample.

The vendor download has no cell-level matrix: it ships a gene x spatial-barcode raw
matrix plus segmentation polygons. This script joins the two.

  1. decode barcode id  SBC:<y_nm>:<x_nm>  ->  (x_um, y_um)      [verified by ..._00_...py]
  2. rasterise the Expanded-5um cell contours into a label mask at RES_UM per pixel
  3. look each barcode up in the mask -> barcode -> cell_id (0 = outside any cell)
  4. stream matrix.mtx.gz and sum barcode counts into (cell, gene)
  5. write <sample>_cells.h5ad  with obs = x, y, area_um2, n_barcodes, nCount, nFeature

Run (one sample at a time; Grade3 is the big one):
  XDG_CACHE_HOME="$TMPDIR/fccache" \
  <CONDA_ROOT>/envs/bin2cell_env/bin/python \
    high_res_ST/novaseq/code/novaseq_01_prep_cells.py --sample IDC_Grade3
"""
import argparse
import gzip
import os
import time

import anndata as ad
import cv2
import numpy as np
import pandas as pd
import scipy.sparse as sp

from novaseq_samples import SAMPLES, PROJ

RES_UM = 0.5          # rasterisation pitch; cells are ~10 um across
CHUNK = 50_000_000    # mtx lines per read
OUT = os.path.join(PROJ, "processed_data", "cells")


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def read_barcode_xy(path):
    """SBC:<y_nm>:<x_nm> -> (x_um, y_um) float32, in the segmentation coordinate frame."""
    df = pd.read_csv(path, sep=":", header=None, usecols=[1, 2],
                     names=["_", "y_nm", "x_nm"], dtype=np.int64)
    return (df.x_nm.values / 1000.0).astype(np.float32), (df.y_nm.values / 1000.0).astype(np.float32)


def rasterise_cells(contour_csv, x0, y0, nx, ny):
    """Label mask (int32, 0 = background) of the cell polygons, plus per-cell pixel area."""
    seg = pd.read_csv(contour_csv, dtype={"cell_id": np.int64,
                                          "vertex_x": np.float32, "vertex_y": np.float32})
    cell_ids = seg.cell_id.values
    px = np.rint((seg.vertex_x.values - x0) / RES_UM).astype(np.int32)
    py = np.rint((seg.vertex_y.values - y0) / RES_UM).astype(np.int32)

    # contour rows are contiguous per cell_id (checked on all three samples)
    starts = np.flatnonzero(np.r_[True, cell_ids[1:] != cell_ids[:-1]])
    ends = np.r_[starts[1:], len(cell_ids)]
    uniq = cell_ids[starts]
    if len(np.unique(uniq)) != len(uniq):
        raise ValueError(f"{contour_csv}: cell_id rows are not contiguous")

    mask = np.zeros((ny, nx), dtype=np.int32)
    for st, en, cid in zip(starts, ends, uniq):
        cv2.fillPoly(mask, [np.stack([px[st:en], py[st:en]], axis=1)], color=int(cid))
    return mask, uniq


def main(sample_name):
    os.makedirs(OUT, exist_ok=True)
    s = next(x for x in SAMPLES if x["name"] == sample_name)
    log("sample:", s["name"])

    bx, by = read_barcode_xy(s["barcodes"])
    log(f"barcodes: {len(bx):,}  x {bx.min():.0f}-{bx.max():.0f}  y {by.min():.0f}-{by.max():.0f}")

    seg_head = pd.read_csv(s["cell_contours"], usecols=["vertex_x", "vertex_y"],
                           dtype=np.float32)
    x0 = min(bx.min(), seg_head.vertex_x.min()) - 1.0
    y0 = min(by.min(), seg_head.vertex_y.min()) - 1.0
    nx = int(np.ceil((max(bx.max(), seg_head.vertex_x.max()) - x0) / RES_UM)) + 2
    ny = int(np.ceil((max(by.max(), seg_head.vertex_y.max()) - y0) / RES_UM)) + 2
    del seg_head
    log(f"mask grid: {ny} x {nx} px @ {RES_UM} um  ({ny * nx * 4 / 1e9:.1f} GB)")

    mask, cell_ids = rasterise_cells(s["cell_contours"], x0, y0, nx, ny)
    log(f"rasterised {len(cell_ids):,} cells")

    # cell_id -> row index; areas straight off the mask
    cid_max = int(cell_ids.max())
    cid2row = np.full(cid_max + 1, -1, dtype=np.int32)
    cid2row[cell_ids] = np.arange(len(cell_ids), dtype=np.int32)
    px_per_cell = np.bincount(mask.ravel(), minlength=cid_max + 1)
    area_um2 = px_per_cell[cell_ids] * RES_UM ** 2

    # barcode -> cell row (-1 = outside every cell)
    ix = np.rint((bx - x0) / RES_UM).astype(np.int32)
    iy = np.rint((by - y0) / RES_UM).astype(np.int32)
    bc_cid = mask[iy, ix]
    del mask, ix, iy, bx, by
    bc_row = np.where(bc_cid > 0, cid2row[np.clip(bc_cid, 0, cid_max)], -1).astype(np.int32)
    in_cell = bc_row >= 0
    log(f"barcodes inside a cell: {in_cell.sum():,} / {len(bc_row):,} "
        f"({100 * in_cell.mean():.1f}%)")
    n_barcodes_per_cell = np.bincount(bc_row[in_cell], minlength=len(cell_ids))

    feats = pd.read_csv(s["features"], sep="\t", header=None,
                        names=["ensembl", "symbol", "type"])
    n_genes = len(feats)

    # ---- stream the mtx, keep only entries whose barcode fell inside a cell ----
    rows, cols, vals = [], [], []
    n_read = 0
    reader = pd.read_csv(s["matrix"], sep=" ", header=None, skiprows=3, chunksize=CHUNK,
                         names=["gene", "bc", "count"], dtype=np.int64, engine="c")
    for k, ch in enumerate(reader):
        g = ch.gene.values - 1
        b = ch.bc.values - 1
        c = ch["count"].values
        r = bc_row[b]
        keep = r >= 0
        rows.append(r[keep].astype(np.int32))
        cols.append(g[keep].astype(np.int32))
        vals.append(c[keep].astype(np.int32))
        n_read += len(ch)
        log(f"  chunk {k}: {n_read:,} entries read, {keep.sum():,} kept")

    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    log(f"building sparse matrix: {len(vals):,} entries -> "
        f"{len(cell_ids):,} cells x {n_genes:,} genes")
    X = sp.coo_matrix((vals, (rows, cols)),
                      shape=(len(cell_ids), n_genes)).tocsr()
    X.sum_duplicates()
    del rows, cols, vals
    log(f"nnz after dedup: {X.nnz:,}")

    # ---- centroids from the polygons (cheap second pass) ----
    seg = pd.read_csv(s["cell_contours"],
                      dtype={"cell_id": np.int64, "vertex_x": np.float32, "vertex_y": np.float32})
    cent = seg.groupby("cell_id")[["vertex_x", "vertex_y"]].mean().loc[cell_ids]
    del seg

    obs = pd.DataFrame(
        dict(cell_id=cell_ids,
             x=cent.vertex_x.values, y=cent.vertex_y.values,
             area_um2=area_um2.astype(np.float32),
             n_barcodes=n_barcodes_per_cell.astype(np.int32),
             nCount_spatial=np.asarray(X.sum(1)).ravel().astype(np.int32),
             nFeature_spatial=np.diff(X.indptr).astype(np.int32)),
        index=[f"{s['short']}_{c}" for c in cell_ids])
    obs["sample"] = s["name"]

    adata = ad.AnnData(X=X, obs=obs, var=feats.set_index("ensembl"))
    adata.obsm["spatial"] = obs[["x", "y"]].values.astype(np.float32)
    out = os.path.join(OUT, f"{s['name']}_cells.h5ad")
    adata.write_h5ad(out, compression="gzip")
    log("DONE ->", out)
    log(f"median nCount {np.median(obs.nCount_spatial):.0f}  "
        f"median nFeature {np.median(obs.nFeature_spatial):.0f}  "
        f"median area {np.median(obs.area_um2):.1f} um2")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, choices=[s["name"] for s in SAMPLES])
    main(ap.parse_args().sample)
