#!/usr/bin/env Rscript
# Build the two side inputs the altered_benign subclone app expects, for the
# IDC_Grade3 ds90k inferCNV run.
#
# The app (atera_brca_analysis/code/altered_benign) reads Xenium-native files:
#   cells_file    parquet with cell_id / x_centroid / y_centroid   (R/io_xenium.R)
#   histology_csv csv with Cell.ID / histology
# NovaSeq has neither, but IDC_Grade3_epi_groups_ds90k.csv already carries x, y and
# cluster_frozen for every cell in the run, so this just reshapes it.
#
# histology is set to the res-0.3 expression cluster ("idc_c0", "idc_c1", ...) rather
# than a pathology call: the app's annotation strip then shows the expression clusters
# next to the CNV tree, and clicking the strip selects one. All idc_* values are
# declared as role "cancer" in the config, so the tree is built on all 90,000 cells.
#
# Env: r4.4 (has arrow).
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/r4.4/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_12_g3_app_inputs.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({ library(arrow) })

sample     <- "IDC_Grade3"
proj       <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir    <- file.path(proj, "processed_data", paste0("infercnv_", sample, "_ds90k"))
groups_csv <- file.path(run_dir, paste0(sample, "_epi_groups_ds90k.csv"))
app_dir    <- file.path(run_dir, "app_inputs")
dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)

g <- read.csv(groups_csv, row.names = 1)

write_parquet(
  data.frame(cell_id = rownames(g), x_centroid = g$x, y_centroid = g$y),
  file.path(app_dir, "cells.parquet"))

hist <- ifelse(g$group == "cancer", paste0("idc_c", g$cluster_frozen), "normal")
write.csv(
  data.frame(Cell.ID = rownames(g), histology = hist),
  file.path(app_dir, "histology_annotation.csv"), row.names = FALSE)

message(sprintf("wrote %d cells -> %s", nrow(g), app_dir))
print(table(hist))
