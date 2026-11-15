#!/usr/bin/env python
"""VisiumHD AP5 — bin2cell nuclei segmentation (2um bins -> single cells).

Step 1/2 of the segmentation chain. Produces the two h5ad files the rest of the
high-resolution analysis is built on:

  <OUT_DIR>/mpp0.5_heprob0.15_gexprob0.2_nms0.5.h5ad            2um bins, carries
      labels_he / labels_he_expanded / labels_gex / labels_joint
  <OUT_DIR>/mpp0.5_heprob0.15_gexprob0.2_nms0.5_countsdata.h5ad cell-level counts
      (bin_to_cell on labels_joint)

Step 2/2 is bin2cell_export_for_seurat.py, which turns these into the CSV/MTX
files read by rctd/run_RCTD_withoutUnassignedCycling.R and rctd/compare_Anno_RCTD.R.

Parameters are the frozen ones (they name the output files): StarDist on the H&E
image at prob_thresh 0.15, StarDist on the gridded gene-expression image at
prob_thresh 0.2 / nms_thresh 0.5, both at mpp 0.5. StarDist is the slow part;
both .npz label files are reused if already present.

Run under the bin2cell env:
  <CONDA_ROOT>/envs/bin2cell_env/bin/python bin2cell_run.py
"""
import os

import bin2cell as b2c
import scanpy as sc

## paths -------------------------------------------------------------------
VISIUM_PATH = "<COLLAB_DATA>/spatial_prostate/oxford/data/HD/H1-JQTF3M2_A1/binned_outputs/square_002um/"
# the image passed to spaceranger --image; the spatial coordinates are based on it
SOURCE_IMAGE = "<COLLAB_DATA>/spatial_prostate/oxford/data/HD/H1-JQTF3M2_A1_HE.tif"
OUT_DIR = "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/segmentation/bin2cell/processed"

## parameters (frozen run) --------------------------------------------------
MPP = 0.5
HE_PROB = 0.15
GEX_PROB = 0.2
NMS = 0.5
STEM = f"mpp{MPP}_heprob{HE_PROB}_gexprob{GEX_PROB}_nms{NMS}"

stardist_dir = os.path.join(OUT_DIR, "stardist")
os.makedirs(stardist_dir, exist_ok=True)

## load 2um bins -----------------------------------------------------------
adata = b2c.read_visium(VISIUM_PATH, source_image_path=SOURCE_IMAGE)
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

## H&E segmentation --------------------------------------------------------
he_img = os.path.join(stardist_dir, "he.tiff")
he_npz = os.path.join(stardist_dir, f"he_prob_{HE_PROB}.npz")
b2c.scaled_he_image(adata, mpp=MPP, save_path=he_img)
b2c.destripe(adata)

if not os.path.exists(he_npz):
    b2c.stardist(image_path=he_img, labels_npz_path=he_npz,
                 stardist_model="2D_versatile_he", prob_thresh=HE_PROB)
b2c.insert_labels(adata, labels_npz_path=he_npz, basis="spatial",
                  spatial_key="spatial_cropped_150_buffer", mpp=MPP,
                  labels_key="labels_he")
b2c.expand_labels(adata, labels_key="labels_he",
                  expanded_labels_key="labels_he_expanded")

## gene-expression segmentation (rescues bins with no nucleus call) --------
gex_img = os.path.join(stardist_dir, "gex.tiff")
gex_npz = os.path.join(stardist_dir, f"gex_nms_{NMS}_prob_{GEX_PROB}.npz")
b2c.grid_image(adata, "n_counts_adjusted", mpp=MPP, sigma=5, save_path=gex_img)

if not os.path.exists(gex_npz):
    b2c.stardist(image_path=gex_img, labels_npz_path=gex_npz,
                 stardist_model="2D_versatile_fluo",
                 prob_thresh=GEX_PROB, nms_thresh=NMS)
b2c.insert_labels(adata, labels_npz_path=gex_npz, basis="array", mpp=MPP,
                  labels_key="labels_gex")

b2c.salvage_secondary_labels(adata, primary_label="labels_he_expanded",
                             secondary_label="labels_gex",
                             labels_key="labels_joint")
adata.write(os.path.join(OUT_DIR, f"{STEM}.h5ad"))

## bins -> cells -----------------------------------------------------------
cdata = b2c.bin_to_cell(adata, labels_key="labels_joint",
                        spatial_keys=["spatial", "spatial_cropped_150_buffer"])
cdata.write(os.path.join(OUT_DIR, f"{STEM}_countsdata.h5ad"))

print(f"DONE bin2cell: {adata.n_obs} bins -> {cdata.n_obs} cells ({STEM})")
