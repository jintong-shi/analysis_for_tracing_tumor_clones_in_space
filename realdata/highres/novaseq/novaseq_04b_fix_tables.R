#!/usr/bin/env Rscript
# Patch for objects produced before the AverageExpression name-mangling fix in
# novaseq_04_cluster.R: AverageExpression() runs make.names() on numeric cluster ids
# ("0" -> "g0"), so the per-cluster n_cells lookup returned NA and the pct_* columns
# with it. Marker scores and compartment calls were unaffected.
#
# Rewrites <sample>_compartment_by_resolution.csv, <sample>_resolution_summary.csv and
# writes <sample>_cell_metadata.csv, straight from the saved object's meta.data.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_04b_fix_tables.R --sample DCIS_IDC_Grade1

suppressWarnings(suppressMessages({ library(Seurat) }))

args   <- commandArgs(trailingOnly = TRUE)
sample <- args[which(args == "--sample") + 1]

PROJ <- "<DATA_ROOT>/high_res_ST/novaseq"
OUT  <- file.path(PROJ, "processed_data", "clustering", sample)
RES  <- seq(0.1, 1.0, by = 0.1)
res_cols <- paste0("spatial_snn_res.", RES)

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

log("reading ", sample, " object")
seu <- readRDS(file.path(OUT, paste0(sample, "_seurat.rds.gz")))
md  <- seu@meta.data

write.csv(cbind(md[, c("x", "y", "nCount_spatial", "nFeature_spatial", "area_um2")],
                md[, res_cols], Embeddings(seu, "umap")),
          file.path(OUT, paste0(sample, "_cell_metadata.csv")))
rm(seu); gc()

asg <- read.csv(file.path(OUT, paste0(sample, "_compartment_by_resolution.csv")))
asg$cluster <- sub("^g", "", asg$cluster)
asg$n_cells <- mapply(function(r, cl)
  sum(md[[paste0("spatial_snn_res.", r)]] == cl), asg$resolution, asg$cluster)
write.csv(asg, file.path(OUT, paste0(sample, "_compartment_by_resolution.csv")), row.names = FALSE)

summ <- do.call(rbind, lapply(split(asg, asg$resolution), function(d) data.frame(
  resolution = d$resolution[1], n_clusters = nrow(d),
  min_margin = min(d$margin), median_margin = median(d$margin),
  n_clusters_margin_lt_0.5 = sum(d$margin < 0.5),
  pct_epithelial = round(100 * sum(d$n_cells[d$compartment == "epithelial"]) / sum(d$n_cells), 1),
  pct_immune     = round(100 * sum(d$n_cells[d$compartment == "immune"]) / sum(d$n_cells), 1),
  pct_stromal    = round(100 * sum(d$n_cells[d$compartment == "stromal"]) / sum(d$n_cells), 1))))
summ <- summ[order(as.numeric(summ$resolution)), ]
write.csv(summ, file.path(OUT, paste0(sample, "_resolution_summary.csv")), row.names = FALSE)
print(summ, row.names = FALSE)
log("DONE -> ", OUT)
