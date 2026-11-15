#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Fig 5 panel d draws two subclone maps (HD epi-meta-cells r40 k=2, regular
# Visium pt10 k=3). prep_r40_k2.rds / prep_visium_k3.rds (subclone_cnv_prep.R)
# hold the CNV matrices and per-unit subclone labels but no coordinates, so the
# panel also reads an id,x,y table per dataset. This script writes those two:
#
#   coords_r40.csv     every meta-cell of the aggregated pseudo-spot object, with
#                      the spatial reduction as x/y. The ids are that object's
#                      colnames (g1, g2, ...) -- the same names the r40 inferCNV
#                      run carries, so they match names(prep_r40_k2.rds$sub).
#   coords_visium.csv  spaceranger full-res pixel coords for every in-tissue
#                      barcode. Panel d needs the benign spots too: they are the
#                      grey background of the Visium map (smap bg_ids).
#
# Do NOT rebuild the meta-cell ids/coords from pseudo_spot/meta_r40_Smax20.csv:
# that mapping CSV and the aggregated *_seu.rds.gz object come from two different
# region-growing runs and cannot be joined by pseudo_id (see the NOTE in
# pseudo_spot/code/plot_metacell_spatial_ncells.R). Joining them silently matches
# ~93% of the ids to the wrong meta-cell.
#
# write.csv (not fwrite) so the output is byte-identical to the frozen inputs.
#
# Run under stcnv-r-shared (SeuratObject):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_subclone_coords.R
# =====================================================================
suppressWarnings(suppressMessages({
  library(data.table); library(SeuratObject)
}))

base  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data"
seu_f <- file.path(base, "pseudo_spot/epi_pseudo_r40_Smax20_seu.rds.gz")
vis_f <- "../visiumHD/histology/visium_tissue_positions.csv"  # spaceranger tissue_positions (AP5_10 / pt10)
OD    <- "../visiumHD/subclone"                               # also holds prep_r40_k2.rds / prep_visium_k3.rds
dir.create(OD, showWarnings = FALSE, recursive = TRUE)

# ids in the prep RDS are authoritative: panel d looks up sub[co$id]. A partial
# match means the two id conventions disagree, which would misplace the clones.
check_ids <- function(ids, rds, what) {
  f <- file.path(OD, rds)
  if (!file.exists(f)) { message(sprintf("[%s] %s not found, skipping id check", what, rds)); return(invisible()) }
  known <- names(readRDS(f)$sub)
  miss  <- setdiff(known, ids)
  if (length(miss) > 0L)
    stop(sprintf("[%s] %d/%d ids in %s have no coords (e.g. %s vs %s) -- id conventions differ",
                 what, length(miss), length(known), rds, miss[1], ids[1]))
  message(sprintf("[%s] all %d %s ids have coords", what, length(known), rds))
}

## ---- HD epithelial meta-cells (r = 40um) ------------------------------------
xy  <- readRDS(seu_f)@reductions$spatial@cell.embeddings   # meta-cells x 2, rownames = g<N>
r40 <- data.frame(id = rownames(xy), x = xy[, 1], y = xy[, 2])
check_ids(r40$id, "prep_r40_k2.rds", "r40")
write.csv(r40, file.path(OD, "coords_r40.csv"), row.names = FALSE)

## ---- regular Visium (pt10) --------------------------------------------------
vp  <- fread(vis_f)
vis <- vp[in_tissue == 1, .(id = barcode, x = pxl_col_in_fullres, y = pxl_row_in_fullres)]
check_ids(vis$id, "prep_visium_k3.rds", "visium")
write.csv(vis, file.path(OD, "coords_visium.csv"), row.names = FALSE)

message(sprintf("DONE -> %s  (coords_r40=%d meta-cells, coords_visium=%d spots)",
                OD, nrow(r40), nrow(vis)))
