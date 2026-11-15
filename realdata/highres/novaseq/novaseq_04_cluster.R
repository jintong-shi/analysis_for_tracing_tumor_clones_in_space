#!/usr/bin/env Rscript
# Deliverable 2a -- Seurat clustering of one NovaSeq X Illumina spatial sample.
#
# LogNormalize -> 2000 variable features -> scale -> PCA (30) -> UMAP (30 PCs) ->
# SNN graph (30 PCs) -> Louvain clusters at resolutions 0.1 .. 1.0.
#
# Then, for every resolution, average marker expression per cluster is z-scored across
# clusters and each cluster is assigned to the epithelial / immune / stromal compartment
# with the highest mean marker z-score. The margin to the runner-up compartment is kept
# so weakly-assigned clusters are visible. Jintong picks the final resolution.
#
# Env: stcnv-r-shared (Seurat 5.1.0 / SeuratObject 5.4.0).
#
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_04_cluster.R --sample IDC_Grade3

suppressWarnings(suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(hdf5r)
}))

args   <- commandArgs(trailingOnly = TRUE)
sample <- args[which(args == "--sample") + 1]

PROJ   <- "<DATA_ROOT>/high_res_ST/novaseq"
CELLS  <- file.path(PROJ, "processed_data", "cells")
OUT    <- file.path(PROJ, "processed_data", "clustering", sample)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RES     <- seq(0.1, 1.0, by = 0.1)
NPC     <- 30
MARKERS <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  immune     = c("CD3D", "CD3E", "JCHAIN", "MZB1", "MS4A2"),
  stromal    = c("COL1A1", "COL1A2", "LUM", "THBS2", "VWF", "PECAM1")
)

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## ---- load ----
log("reading ", sample)
counts <- Read10X_h5(file.path(CELLS, paste0(sample, "_cells_filtered.h5")))
obs <- read.csv(file.path(CELLS, paste0(sample, "_cells_filtered_obs.csv")), row.names = 1)
obs <- obs[colnames(counts), ]
seu <- CreateSeuratObject(counts, project = sample, meta.data = obs, assay = "spatial")
log("cells: ", ncol(seu), "  genes: ", nrow(seu))

## ---- standard pipeline ----
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
log("PCA")
seu <- RunPCA(seu, npcs = NPC, verbose = FALSE)
log("UMAP")
seu <- RunUMAP(seu, dims = 1:NPC, verbose = FALSE)
log("SNN graph")
seu <- FindNeighbors(seu, dims = 1:NPC, verbose = FALSE)
for (r in RES) {
  log("clustering res ", r)
  seu <- FindClusters(seu, resolution = r, verbose = FALSE)
}
res_cols <- paste0("spatial_snn_res.", RES)

saveRDS(seu, file.path(OUT, paste0(sample, "_seurat.rds.gz")), compress = "gzip")
write.csv(cbind(seu@meta.data[, c("x", "y", "nCount_spatial", "nFeature_spatial", "area_um2")],
                seu@meta.data[, res_cols], Embeddings(seu, "umap")),
          file.path(OUT, paste0(sample, "_cell_metadata.csv")))
log("object saved")

## ---- compartment assignment per resolution ----
genes <- unlist(MARKERS, use.names = FALSE)
genes <- genes[genes %in% rownames(seu)]
missing <- setdiff(unlist(MARKERS, use.names = FALSE), genes)
if (length(missing)) log("markers absent from matrix: ", paste(missing, collapse = ", "))

assign_tbl <- list(); expr_tbl <- list()
for (rc in res_cols) {
  Idents(seu) <- rc
  avg <- AverageExpression(seu, features = genes, assays = "spatial",
                           layer = "data", verbose = FALSE)$spatial
  # AverageExpression runs make.names() on the idents, so numeric cluster ids come
  # back as "g0", "g1", ...; columns are in level order, so restore the real names.
  colnames(avg) <- levels(seu[[rc]][, 1])
  avg <- log1p(avg)
  z <- t(scale(t(avg)))                       # z-score each gene across clusters
  z[is.na(z)] <- 0
  score <- sapply(MARKERS, function(gs) colMeans(z[intersect(gs, rownames(z)), , drop = FALSE]))
  score <- matrix(score, nrow = ncol(z), dimnames = list(colnames(z), names(MARKERS)))
  best  <- colnames(score)[max.col(score)]
  srt   <- t(apply(score, 1, sort, decreasing = TRUE))
  assign_tbl[[rc]] <- data.frame(
    resolution = sub(".*res\\.", "", rc), cluster = rownames(score),
    n_cells = as.integer(table(seu[[rc]][, 1])[rownames(score)]),
    compartment = best, margin = round(srt[, 1] - srt[, 2], 3),
    round(score, 3), row.names = NULL)
  expr_tbl[[rc]] <- data.frame(resolution = sub(".*res\\.", "", rc),
                               cluster = colnames(avg), t(round(avg, 3)), row.names = NULL)
}
assign_df <- do.call(rbind, assign_tbl)
write.csv(assign_df, file.path(OUT, paste0(sample, "_compartment_by_resolution.csv")), row.names = FALSE)
write.csv(do.call(rbind, expr_tbl), file.path(OUT, paste0(sample, "_marker_avgexpr.csv")), row.names = FALSE)

## ---- summary per resolution: n clusters, weakest margin, compartment split ----
summ <- do.call(rbind, lapply(split(assign_df, assign_df$resolution), function(d) data.frame(
  resolution = d$resolution[1], n_clusters = nrow(d),
  min_margin = min(d$margin), median_margin = median(d$margin),
  n_clusters_margin_lt_0.5 = sum(d$margin < 0.5),
  pct_epithelial = round(100 * sum(d$n_cells[d$compartment == "epithelial"]) / sum(d$n_cells), 1),
  pct_immune     = round(100 * sum(d$n_cells[d$compartment == "immune"]) / sum(d$n_cells), 1),
  pct_stromal    = round(100 * sum(d$n_cells[d$compartment == "stromal"]) / sum(d$n_cells), 1))))
write.csv(summ, file.path(OUT, paste0(sample, "_resolution_summary.csv")), row.names = FALSE)
print(summ, row.names = FALSE)

## ---- figures ----
save_both <- function(p, base, w, h) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  log("DONE -> ", basename(base), ".{png,pdf}")
}

md <- seu@meta.data
md$umap1 <- Embeddings(seu, "umap")[, 1]; md$umap2 <- Embeddings(seu, "umap")[, 2]

umap_grid <- wrap_plots(lapply(res_cols, function(rc) {
  ggplot(md, aes(umap1, umap2, colour = .data[[rc]])) +
    geom_point(size = 0.1, stroke = 0) + theme_void() +
    ggtitle(paste0("res ", sub(".*res\\.", "", rc), " (", nlevels(md[[rc]]), " clusters)")) +
    theme(legend.position = "none", plot.title = element_text(size = 8))
}), ncol = 5)
save_both(umap_grid, file.path(OUT, paste0(sample, "_umap_by_resolution")), 16, 7)

spatial_grid <- wrap_plots(lapply(res_cols, function(rc) {
  ggplot(md, aes(x, y, colour = .data[[rc]])) +
    geom_point(size = 0.1, stroke = 0) + coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(paste0("res ", sub(".*res\\.", "", rc))) +
    theme(legend.position = "none", plot.title = element_text(size = 8))
}), ncol = 5)
save_both(spatial_grid, file.path(OUT, paste0(sample, "_spatial_by_resolution")), 16, 9)

for (rc in res_cols) {
  Idents(seu) <- rc
  p <- DotPlot(seu, features = genes, cluster.idents = TRUE) +
    RotatedAxis() + ggtitle(paste0(sample, " -- ", rc)) +
    theme(axis.text = element_text(size = 7), plot.title = element_text(size = 9))
  save_both(p, file.path(OUT, paste0(sample, "_dotplot_", rc)), 10, 2 + 0.18 * nlevels(Idents(seu)))
}

log("ALL DONE -> ", OUT)
