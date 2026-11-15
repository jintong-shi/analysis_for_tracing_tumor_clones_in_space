#!/usr/bin/env Rscript
# benign_ref_strictqc_infcnv.R
#
# Run 4 - rerun of the run2 benign-reference inferCNV, but with a stricter QC
# cutoff applied to the epithelial cells.
#
# Differs from run2/4.qc_benign_ref_infcnv.R only in that:
#   - it starts from the base epithelial object (xenium_epi.rds.gz, nCount >= 100)
#     and applies the stricter QC (nCount > 1000 & nFeature > 800) here, and
#   - the infercnv::run() call is enabled (in script 4 it is commented out).
#
# Everything else mirrors script 4:
#   1. Load epithelial object, apply QC subset.
#   2. Add histology group (normal_1/normal_2 = normal, else cancer).
#   3. Benign (histology "normal") reference inferCNV on all QC-passing epi cells.
#   4. Reference/observation ComplexHeatmap.
#
# inferCNV parameters (identical to script 4 / run3):
#   cutoff = 0.1, denoise = TRUE, HMM = FALSE, analysis_mode = "cells",
#   cluster_by_groups = FALSE, write_phylo = TRUE, num_threads = 4,
#   chr_exclude = "chrM", gene_order_file = gencode_v38.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(infercnv)
  library(readr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)

## ---- config -----------------------------------------------------
epi_path      <- "<DATA_ROOT>/high_res_ST/atera/processed_data/clustering/xenium_epi.rds.gz"
histology_csv <- "<DATA_ROOT>/high_res_ST/atera/data/histology/histology_annotation.csv"
output_dir    <- "<DATA_ROOT>/high_res_ST/atera/processed_data/run4_benign_ref_strictqc"
gene_pos_file <- "<REF_ROOT>/gencode_v38_gene_pos.txt"

assay         <- "Xenium"
counts_layer  <- "counts"

# stricter QC cutoff (this is the only analysis change vs run2)
min_count     <- 1000
min_feature   <- 800

cutoff        <- 0.1
ncores        <- 4

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- reference / observation ComplexHeatmap ---------------------
## (verbatim from run2/4.qc_benign_ref_infcnv.R)
plot_infercnv_ref_obs_heatmap <- function(
    infercnv_obj,
    output_file,
    max_cells_per_panel = 5000,
    seed = 42,
    cnv_bounds = c(0.853, 1.0, 1.157),
    cnv_colors = c("#2166AC", "white", "#B2182B"),
    width = 3000,
    height = 1800,
    res = 150,
    raster_quality = 2
) {
  if (is.character(infercnv_obj)) {
    infercnv_obj <- readRDS(infercnv_obj)
  }

  cnv_mat <- infercnv_obj@expr.data

  get_cell_indices <- function(grouped_cell_indices, panel_name) {
    cell_indices <- unlist(grouped_cell_indices, use.names = FALSE)
    if (length(cell_indices) == 0L) {
      stop("No ", panel_name, " cell indices found in the infercnv object.")
    }

    if (is.character(cell_indices)) {
      cell_indices <- match(cell_indices, colnames(cnv_mat))
    }

    cell_indices <- unique(as.integer(cell_indices))
    cell_indices <- cell_indices[!is.na(cell_indices)]
    cell_indices <- cell_indices[cell_indices >= 1L & cell_indices <= ncol(cnv_mat)]
    if (length(cell_indices) == 0L) {
      stop("No valid ", panel_name, " cell indices matched expr.data columns.")
    }

    cell_indices
  }

  sample_cells <- function(cell_indices) {
    if (!is.null(max_cells_per_panel) && length(cell_indices) > max_cells_per_panel) {
      cell_indices <- sample(cell_indices, max_cells_per_panel)
    }
    cell_indices
  }

  ref_indices <- get_cell_indices(infercnv_obj@reference_grouped_cell_indices, "reference")
  obs_indices <- get_cell_indices(infercnv_obj@observation_grouped_cell_indices, "observation")

  if (!is.null(seed)) {
    set.seed(seed)
  }
  ref_keep <- sample_cells(ref_indices)
  obs_keep <- sample_cells(obs_indices)

  col_fun <- colorRamp2(cnv_bounds, cnv_colors)
  heatmap_args <- list(
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_quality = raster_quality
  )

  ref_heatmap <- do.call(
    Heatmap,
    c(
      list(
        matrix = t(cnv_mat[, ref_keep, drop = FALSE]),
        name = "CNV",
        column_title = paste0("Reference cells (n = ", length(ref_keep), ")")
      ),
      heatmap_args
    )
  )
  obs_heatmap <- do.call(
    Heatmap,
    c(
      list(
        matrix = t(cnv_mat[, obs_keep, drop = FALSE]),
        name = "CNV_observation",
        column_title = paste0("Observation cells (n = ", length(obs_keep), ")"),
        show_heatmap_legend = FALSE
      ),
      heatmap_args
    )
  )

  png(output_file, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  draw(ref_heatmap %v% obs_heatmap, heatmap_legend_side = "right")

  invisible(list(
    reference_cells = colnames(cnv_mat)[ref_keep],
    observation_cells = colnames(cnv_mat)[obs_keep],
    output_file = output_file
  ))
}

## ---- 1. load epithelial object + apply stricter QC --------------
message("Loading epithelial object: ", epi_path)
xenium.epi <- readRDS(epi_path)
xenium.epi <- UpdateSeuratObject(xenium.epi)
n_pre <- ncol(xenium.epi)

xenium.epi <- subset(xenium.epi,
                     subset = nCount_Xenium > min_count & nFeature_Xenium > min_feature)
message(sprintf("QC nCount > %d & nFeature > %d: kept %d / %d epithelial cells",
                min_count, min_feature, ncol(xenium.epi), n_pre))

## ---- 2. add histology group (normal / cancer) -------------------
histology <- read.csv(histology_csv)
histology <- tibble::column_to_rownames(histology, "Cell.ID")
histology <- histology %>%
  dplyr::mutate(group = ifelse(histology %in% c("normal_1", "normal_2"), "normal", "cancer"))
xenium.epi <- AddMetaData(xenium.epi, dplyr::select(histology, group))

n_na <- sum(is.na(xenium.epi$group))
if (n_na > 0) {
  message(sprintf("WARNING: %d QC-passing epi cells have no histology annotation (group = NA).", n_na))
}
message("group table (QC-passing epi cells):")
print(table(xenium.epi$group, useNA = "ifany"))

saveRDS(xenium.epi,
        file.path(output_dir, "xenium_with_anno_strictqc.rds.gz"),
        compress = "gzip")

## ---- 3. benign-reference inferCNV on all QC-passing epi ---------
infcnv_dir <- file.path(output_dir, "infcnv")
dir.create(infcnv_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(infcnv_dir, "run.final.infercnv_obj"))) {
  group.dat <- data.frame(row.names = rownames(xenium.epi@meta.data),
                          group = xenium.epi$group)
  expr_mat  <- GetAssayData(xenium.epi, assay = assay, layer = counts_layer)

  message("Running benign-reference inferCNV (ref = normal) ...")
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = expr_mat,
    gene_order_file   = gene_pos_file,
    annotations_file  = group.dat,
    delim             = "\t",
    ref_group_names   = "normal",
    chr_exclude       = "chrM"
  )
  infercnv::run(
    obj, cutoff = cutoff, out_dir = infcnv_dir,
    denoise = TRUE, HMM = FALSE,
    cluster_by_groups = FALSE,
    write_phylo = TRUE,
    analysis_mode = "cells",
    num_threads = ncores,
    no_plot = TRUE
  )
}
final_obj <- file.path(infcnv_dir, "run.final.infercnv_obj")
message("inferCNV done: ", final_obj)

## ---- 4. reference / observation heatmaps ------------------------
plot_infercnv_ref_obs_heatmap(
  infercnv_obj = final_obj,
  output_file  = file.path(infcnv_dir, "infercnv_complexheatmap_downsample5000.png"),
  max_cells_per_panel = 5000
)
plot_infercnv_ref_obs_heatmap(
  infercnv_obj = final_obj,
  output_file  = file.path(infcnv_dir, "infercnv_complexheatmap.png"),
  max_cells_per_panel = NULL,
  height = 2200
)
message("Heatmaps written to ", infcnv_dir)
message("Done.")
