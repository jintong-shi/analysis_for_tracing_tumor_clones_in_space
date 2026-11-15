#!/usr/bin/env Rscript
# novaseq_09_g3_infcnv_heatmap.R
#
# Reference / observation ComplexHeatmap for the IDC_Grade3 benign-reference inferCNV
# run (novaseq_08_g3_benign_ref_infcnv.R). Split out because no conda env carries both
# infercnv and ComplexHeatmap.
#
# The plotting function is verbatim from
# opensource/benchmark_spatialCNV/realdata/highres/atera/benign_ref_strictqc_infcnv.R;
# the second call is capped at 50,000 cells per panel instead of NULL, because
# clustering 520k rows in ComplexHeatmap is not feasible.
#
# Env: stcnv-r-tuscan.
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_09_g3_infcnv_heatmap.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)

sample     <- "IDC_Grade3"
proj       <- "<DATA_ROOT>/high_res_ST/novaseq"
# optional first argument: the infcnv output dir to plot (e.g. the ds90k run)
args       <- commandArgs(trailingOnly = TRUE)
infcnv_dir <- if (length(args) > 0) args[1] else
  file.path(proj, "processed_data", paste0("infercnv_", sample), "infcnv")
final_obj  <- file.path(infcnv_dir, "run.final.infercnv_obj")

## ---- reference / observation ComplexHeatmap ---------------------
## (verbatim from atera/benign_ref_strictqc_infcnv.R)
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


if (!file.exists(final_obj)) {
  stop("inferCNV has not finished yet: ", final_obj, " is missing.")
}

plot_infercnv_ref_obs_heatmap(
  infercnv_obj = final_obj,
  output_file  = file.path(infcnv_dir, "infercnv_complexheatmap_downsample5000.png"),
  max_cells_per_panel = 5000
)
plot_infercnv_ref_obs_heatmap(
  infercnv_obj = final_obj,
  output_file  = file.path(infcnv_dir, "infercnv_complexheatmap_downsample50000.png"),
  max_cells_per_panel = 50000,
  height = 2200
)
message("Heatmaps written to ", infcnv_dir)
