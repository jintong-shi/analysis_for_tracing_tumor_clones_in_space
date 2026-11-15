#!/usr/bin/env python
"""Export each cell x gene h5ad as a 10x-format .h5 that Seurat's Read10X_h5 can load.

A text .mtx would be ~25 GB for Grade3 (1.2e9 nonzeros); the 10x HDF5 layout is the
compact equivalent. Our AnnData X is cells x genes CSR, which is bit-for-bit the CSC
of the genes x cells matrix 10x stores, so no transpose is needed.

Cells are filtered to the deliverable-1 thresholds and genes to those detected in at
least MIN_CELLS cells. The unfiltered h5ad stays the source of truth (inferCNV will
want all genes back).

Run:
  <CONDA_ROOT>/envs/bin2cell_env/bin/python \
    high_res_ST/novaseq/code/novaseq_03_export_10x.py
"""
import os

import anndata as ad
import h5py
import numpy as np

from novaseq_samples import SAMPLES, PROJ

CELLS = os.path.join(PROJ, "processed_data", "cells")
MIN_COUNT, MIN_FEATURE, MIN_CELLS = 50, 25, 10


def write_10x_h5(path, X, barcodes, gene_ids, gene_names):
    """X: cells x genes CSR -> /matrix stored as genes x cells CSC."""
    if os.path.exists(path):
        os.remove(path)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        g.create_dataset("data", data=X.data.astype(np.int32), compression="gzip", compression_opts=1)
        g.create_dataset("indices", data=X.indices.astype(np.int64), compression="gzip", compression_opts=1)
        g.create_dataset("indptr", data=X.indptr.astype(np.int64), compression="gzip", compression_opts=1)
        g.create_dataset("shape", data=np.array([X.shape[1], X.shape[0]], dtype=np.int32))
        fe = g.create_group("features")
        fe.create_dataset("id", data=np.array(gene_ids, dtype="S"))
        fe.create_dataset("name", data=np.array(gene_names, dtype="S"))
        fe.create_dataset("feature_type", data=np.array(["Gene Expression"] * X.shape[1], dtype="S"))
        fe.create_dataset("genome", data=np.array(["GRCh38"] * X.shape[1], dtype="S"))
        fe.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


def main():
    for s in SAMPLES:
        src = os.path.join(CELLS, f"{s['name']}_cells.h5ad")
        a = ad.read_h5ad(src)
        keep_c = (a.obs.nCount_spatial >= MIN_COUNT) & (a.obs.nFeature_spatial >= MIN_FEATURE)
        a = a[keep_c.values].copy()
        n_cells_per_gene = np.diff(a.X.tocsc().indptr)
        keep_g = n_cells_per_gene >= MIN_CELLS
        a = a[:, keep_g].copy()
        print(f"{s['name']}: {a.n_obs:,} cells x {a.n_vars:,} genes, {a.X.nnz:,} nnz", flush=True)

        # gene symbols are not unique in this reference; Seurat needs unique rownames
        names = a.var["symbol"].astype(str).values
        write_10x_h5(os.path.join(CELLS, f"{s['name']}_cells_filtered.h5"),
                     a.X.tocsr(), a.obs_names.tolist(), a.var_names.tolist(), names)
        a.obs.to_csv(os.path.join(CELLS, f"{s['name']}_cells_filtered_obs.csv"))
        print("  DONE ->", f"{s['name']}_cells_filtered.h5", flush=True)


if __name__ == "__main__":
    main()
