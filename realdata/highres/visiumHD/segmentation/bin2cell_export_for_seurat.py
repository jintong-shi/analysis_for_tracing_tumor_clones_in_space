#!/usr/bin/env python
"""VisiumHD AP5 — export the bin2cell h5ad files for the R side.

Step 2/2 of the segmentation chain (step 1 is bin2cell_run.py). Writes the flat
files the downstream R scripts read; nothing is recomputed here.

From the 2um bin object:
  <SEG_DIR>/metadata_labels_joint.csv   obs of every 2um bin, incl. labels_joint.
      Read by rctd/compare_Anno_RCTD.R.

From the cell-level object:
  <SEG_DIR>/{matrix.mtx,features.tsv,barcodes.tsv}  counts, genes x cells
  <SEG_DIR>/metadata.csv                            per-cell obs
  <SEG_DIR>/spatial.csv                             per-cell x/y
      These five are assembled into bin2cell_seu_counts.rds.gz by the
      CreateSeuratObject block at the top of rctd/run_RCTD_withoutUnassignedCycling.R.

Run under the bin2cell env:
  <CONDA_ROOT>/envs/bin2cell_env/bin/python bin2cell_export_for_seurat.py
"""
import os

import anndata as ad
import pandas as pd
from scipy.io import mmwrite

SEG_DIR = "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/segmentation"
PROC_DIR = os.path.join(SEG_DIR, "bin2cell", "processed")
STEM = "mpp0.5_heprob0.15_gexprob0.2_nms0.5"

CELL_META_COLS = ["object_id", "bin_count", "array_row", "array_col",
                  "labels_joint_source"]

## 2um bins -> metadata_labels_joint.csv -----------------------------------
adata = ad.read_h5ad(os.path.join(PROC_DIR, f"{STEM}.h5ad"))
adata.obs.to_csv(os.path.join(SEG_DIR, "metadata_labels_joint.csv"))
print(f"metadata_labels_joint.csv: {adata.n_obs} bins")
del adata

## cells -> mtx + metadata + coordinates -----------------------------------
cdata = ad.read_h5ad(os.path.join(PROC_DIR, f"{STEM}_countsdata.h5ad"))

mmwrite(os.path.join(SEG_DIR, "matrix.mtx"), cdata.X.T.tocoo())
pd.Series(cdata.var_names).to_csv(os.path.join(SEG_DIR, "features.tsv"),
                                  sep="\t", index=False, header=False)
pd.Series(cdata.obs_names).to_csv(os.path.join(SEG_DIR, "barcodes.tsv"),
                                  sep="\t", index=False, header=False)

cdata.obs[CELL_META_COLS].to_csv(os.path.join(SEG_DIR, "metadata.csv"))
pd.DataFrame(cdata.obsm["spatial"], index=cdata.obs_names,
             columns=["x", "y"]).to_csv(os.path.join(SEG_DIR, "spatial.csv"))

print(f"cell-level export: {cdata.n_obs} cells x {cdata.n_vars} genes")