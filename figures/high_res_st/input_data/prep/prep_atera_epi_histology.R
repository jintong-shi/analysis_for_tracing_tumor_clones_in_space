#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Fig 5 panel e ("atera epi spatial") colours the Xenium epithelial cells by
# histology, Cancer (DCIS) vs Normal. Both frozen inputs come from the same run4
# benign-ref inferCNV object that panels g/h use, so panel e is guaranteed to
# show exactly the cell set the subclone panels were computed on:
#
#   cells_coords.csv              cell_id, x_centroid, y_centroid (all Xenium cells,
#                                 the panel joins the labelled subset out of it)
#   epi_cancer_normal_labels.csv  cell_id, label
#
# The label IS the run4 reference/observation split: that run was set up with
# ref_group_names = "normal" over the histology group
# (normal_1/normal_2 -> "normal", else "cancer"), see
# realdata/highres/atera/benign_ref_strictqc_infcnv.R lines 172 and 201.
# So reference cells are the histology-normal epithelium and observation cells
# are the histology-cancer (DCIS) epithelium -- no re-derivation needed here.
#
# Object slots are read via attr() (same trick as prepare_atera_subclones.R) so
# this does not need the infercnv package.
#
# epi_cancer_normal_labels.csv reproduces the frozen copy byte for byte.
# cells_coords.csv is value-identical but not byte-identical: the frozen copy was
# written from pandas, whose shortest round-trip float repr keeps 17 significant
# digits, while R writes 15. The Xenium centroids are float32 (~7 digits), so the
# difference is below the precision of the source data.
#
# NOT run end-to-end here (no data in this code-only checkout); parse-checked only.
#
# Run:
#   Rscript prep_atera_epi_histology.R <infercnv_obj> <cells_parquet> [out_dir]
#   <infercnv_obj>  : run4 run.final.infercnv_obj (benign-ref, epithelial)
#   <cells_parquet> : Xenium data/cells.parquet (cell_id, x_centroid, y_centroid)
#   <out_dir>       : e.g. ../atera  (figures/high_res_st/input_data/atera)
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow)
})

args       <- commandArgs(trailingOnly = TRUE)
obj_path   <- args[1]
cells_file <- args[2]
out_dir    <- if (length(args) >= 3) args[3] else "../atera"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- run4 observation / reference cells -------------------------------------
obj <- suppressWarnings(readRDS(obj_path))
ed  <- attr(obj, "expr.data")
obs <- attr(obj, "observation_grouped_cell_indices")
ref <- attr(obj, "reference_grouped_cell_indices")
if (is.null(ed) || is.null(obs) || is.null(ref))
  stop("unexpected inferCNV object layout: ", obj_path)

cells     <- colnames(ed)
obs_cells <- cells[unlist(obs, use.names = FALSE)]
ref_cells <- cells[unlist(ref, use.names = FALSE)]

lab <- data.table(
  cell_id = c(obs_cells, ref_cells),
  label   = rep(c("Cancer (DCIS)", "Normal"), c(length(obs_cells), length(ref_cells))))
message("label breakdown:"); print(table(lab$label))

## ---- per-cell centroids (whole section, panel e inner-joins on cell_id) -----
cc <- as.data.table(read_parquet(
  cells_file, col_select = c("cell_id", "x_centroid", "y_centroid")))

# the join is an inner one, so a coords gap would silently drop labelled cells
missing <- setdiff(lab$cell_id, cc$cell_id)
if (length(missing) > 0L)
  stop(sprintf("%d labelled cells have no centroid in %s (e.g. %s)",
               length(missing), basename(cells_file), missing[1]))

fwrite(cc, file.path(out_dir, "cells_coords.csv"))
write.csv(lab, file.path(out_dir, "epi_cancer_normal_labels.csv"), row.names = FALSE)

message("DONE -> ", out_dir,
        sprintf("  (cancer=%d, normal=%d)", length(obs_cells), length(ref_cells)))