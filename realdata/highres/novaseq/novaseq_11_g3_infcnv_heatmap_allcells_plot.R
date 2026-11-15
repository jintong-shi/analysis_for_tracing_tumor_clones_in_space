#!/usr/bin/env Rscript
# All-cell IDC_Grade3 ds90k inferCNV heatmap - render (step 2/2).
#
# Stacked ComplexHeatmap: benign reference on top (n=476, unbinned), then all 90,000
# cancer cells below as bin-averaged rows, split by cluster_frozen. No row clustering
# anywhere. Columns = genes genome-ordered, split by chromosome.
# Colour scale matches novaseq_09_g3_infcnv_heatmap.R (0.853 / 1.0 / 1.157).
#
# Env: stcnv-r-tuscan (only env with ComplexHeatmap), after the prep step.
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_11_g3_infcnv_heatmap_allcells_plot.R
# ------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
}))

## --kmeans <k>: render the RDS the prep step wrote with --labels (CNV k-means).
args    <- commandArgs(trailingOnly = TRUE)
km      <- if ("--kmeans" %in% args) args[[which(args == "--kmeans") + 1L]] else NA_character_
tag     <- if ("--tag" %in% args) args[[which(args == "--tag") + 1L]] else
           if (is.na(km)) "allcells" else sprintf("allcells_kmeans_k%s", km)
grp_by  <- if (is.na(km) && !("--tag" %in% args)) "cluster_frozen" else "CNV cluster"

sample  <- "IDC_Grade3"
proj    <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir <- file.path(proj, "processed_data", paste0("infercnv_", sample, "_ds90k"))
in_rds  <- file.path(run_dir, "outputs", paste0(sample, "_ds90k_", tag, "_heatmap_input.rds"))
out_png <- file.path(run_dir, "outputs", paste0(sample, "_ds90k_", tag, "_heatmap.png"))

t0   <- Sys.time()
tlog <- function(m) message(sprintf("[t+%6.1fs] %s",
                    as.numeric(difftime(Sys.time(), t0, units = "secs")), m))
d <- readRDS(in_rds)
tlog(sprintf("loaded RDS (ref=%d, obs=%d -> %d binned rows, genes=%d)",
             d$n_ref, d$n_obs, nrow(d$mat_bin), ncol(d$mat_bin)))

col_fun  <- colorRamp2(c(0.853, 1.0, 1.157), c("#2166AC", "white", "#B2182B"))
grp_lv   <- c("Benign", d$lab_levels)
grp_cols <- setNames(c("#55A868", hcl.colors(length(d$lab_levels), "Dark 3")), grp_lv)

common <- function(...) list(
  col = col_fun, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  column_split = d$chr_genes, column_gap = unit(0.5, "mm"),
  column_title_gp = gpar(fontsize = 8), border = TRUE,
  use_raster = TRUE, raster_quality = 2,
  row_title_rot = 0, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  ...
)

ht_ref <- do.call(Heatmap, c(list(
  d$mat_ref, name = "inferCNV",
  row_title = sprintf("Benign\n(ref, n=%d)", d$n_ref),
  cluster_rows = FALSE, height = unit(1, "null"),
  left_annotation = rowAnnotation(
    Group = factor(rep("Benign", d$n_ref), levels = grp_lv),
    col = list(Group = grp_cols), show_annotation_name = FALSE),
  heatmap_legend_param = list(title = "inferCNV", at = c(0.853, 1.0, 1.157))
), common()))

row_titles <- sprintf("%s\n(n=%d)", d$lab_levels, d$n_per_lab[d$lab_levels])
ht_obs <- do.call(Heatmap, c(list(
  d$mat_bin, name = "inferCNV_obs", show_heatmap_legend = FALSE,
  row_split = factor(d$bin_lab, levels = d$lab_levels),
  row_title = row_titles,
  cluster_rows = FALSE, cluster_row_slices = FALSE, height = unit(6, "null"),
  left_annotation = rowAnnotation(
    Group = factor(d$bin_lab, levels = grp_lv),
    col = list(Group = grp_cols), show_annotation_name = FALSE)
), common()))

tlog("heatmap objects built; drawing to raster PNG ...")
png(out_png, width = 3000, height = 2200, res = 150)
draw(ht_ref %v% ht_obs, ht_gap = unit(4, "mm"),
     column_title = sprintf(
       "%s ds90k inferCNV - benign ref (n=%d) + all cancer cells by %s (n=%d)\n(cancer rows = mean of %d cells for rendering)",
       sample, d$n_ref, grp_by, d$n_obs, d$bin),
     column_title_gp = gpar(fontsize = 16, fontface = "bold"),
     merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
tlog("draw complete")
message("DONE ", out_png)
