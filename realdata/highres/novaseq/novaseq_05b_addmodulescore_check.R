#!/usr/bin/env Rscript
# One-off check: does the compartment call survive a standard AddModuleScore?
#
# novaseq_04_cluster.R assigns each cluster to epithelial / immune / stromal from
# per-cluster marker z-scores: AverageExpression on the data layer, log1p again,
# z-score each gene across clusters, mean per module, argmax. Two things in that
# recipe are worth testing:
#   - the second log1p is applied to values that are already log-normalized, which
#     compresses the high-expressing genes (the KRTs) more than the low ones;
#   - z-scoring across 12 clusters makes every score relative to this section, so a
#     cluster can "win" a module while expressing none of it (cluster 9 is called
#     epi with all three scores negative).
#
# Seurat's AddModuleScore does it per cell instead, against expression-matched
# control genes, so it is absolute rather than relative. This script runs it with
# the same three gene sets and reports where the two disagree. It changes nothing;
# the pipeline and the figures are untouched.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05b_addmodulescore_check.R

suppressWarnings(suppressMessages({ library(Seurat); library(data.table); library(ggplot2) }))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# the same panel as novaseq_04_cluster.R:34-36
MARKERS <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  immune     = c("CD3D", "CD3E", "JCHAIN", "MZB1", "MS4A2"),
  stromal    = c("COL1A1", "COL1A2", "LUM", "THBS2", "VWF", "PECAM1"))
AMBIG_MARGIN <- 0.5          # novaseq_05_annotate.R:31
MAP <- c(epithelial = "epi", immune = "imm", stromal = "stroma")

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

log("reading the annotated object (7 GB gz, a few minutes)")
seu <- readRDS(file.path(CLUST, paste0(SAMPLE, "_seurat_annotated_res", RES, ".rds.gz")))
DefaultAssay(seu) <- "spatial"
log("cells: ", ncol(seu), "  genes: ", nrow(seu))

present <- lapply(MARKERS, intersect, rownames(seu))
stopifnot(lengths(present) == lengths(MARKERS))

# Seurat 5.1.0's AddModuleScore calls GetAssayData(slot = ), which SeuratObject
# 5.4.0 made defunct, so the shipped function errors in this env. Rather than touch
# a shared conda env, its body is reproduced here verbatim from
# `print(Seurat::AddModuleScore)` with the one substitution layer = instead of
# slot = : same pool, same binning, same control sampling order, same seed, so the
# RNG draws and therefore the scores are the ones Seurat would have produced.
add_module_score <- function(object, features, nbin = 24, ctrl = 100, seed = 1) {
  set.seed(seed)
  assay.data <- GetAssayData(object, assay = "spatial", layer = "data")
  pool <- rownames(object)
  data.avg <- Matrix::rowMeans(assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(data.avg + rnorm(length(data.avg)) / 1e30,
                                  n = nbin, labels = FALSE, right = FALSE)
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector("list", length(features))
  for (i in seq_along(features)) {
    for (g in features[[i]])
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(sample(data.cut[which(data.cut == data.cut[g])],
                                      size = ctrl, replace = FALSE)))
  }
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- t(sapply(ctrl.use, function(g) Matrix::colMeans(assay.data[g, ])))
  feat.scores <- t(sapply(features, function(g) Matrix::colMeans(assay.data[g, , drop = FALSE])))
  out <- as.data.frame(t(feat.scores - ctrl.scores))
  colnames(out) <- names(features); rownames(out) <- colnames(object)
  out
}

log("AddModuleScore (Seurat 5.1.0 algorithm, nbin 24, ctrl 100, seed 1)")
sc <- add_module_score(seu, present)
md <- as.data.table(cbind(cell = rownames(seu@meta.data),
                          seu@meta.data[, c("cluster_frozen", "major_celltype")], sc))
md[, cluster := as.character(cluster_frozen)]

# per-cluster mean of the per-cell scores, then the same argmax + margin rule
new <- md[, c(.(n_cells = .N, frozen_call = as.character(major_celltype[1])),
              lapply(.SD, mean)), by = cluster, .SDcols = names(MARKERS)]
S <- as.matrix(new[, names(MARKERS), with = FALSE])
new[, ams_compartment := colnames(S)[max.col(S)]]
new[, ams_call := unname(MAP[ams_compartment])]
srt <- t(apply(S, 1, sort, decreasing = TRUE))
new[, ams_margin := round(srt[, 1] - srt[, 2], 3)]
new[, ams_ambiguous := ams_margin < AMBIG_MARGIN]

old <- fread(file.path(CLUST, paste0(SAMPLE, "_annotation_res", RES, ".csv")))
old[, cluster := as.character(cluster)]
cmp <- merge(new, old[, .(cluster, old_call = major_celltype, old_margin = margin,
                          old_ambiguous = ambiguous)], by = "cluster")
cmp[, agrees := ams_call == old_call]
setcolorder(cmp, c("cluster", "n_cells", "old_call", "ams_call", "agrees",
                   "old_margin", "ams_margin", "old_ambiguous", "ams_ambiguous"))
cmp <- cmp[order(-n_cells)]
fwrite(cmp, file.path(OUT, sprintf("%s_res%s_addmodulescore_vs_frozen.csv", SAMPLE, RES)))
fwrite(md[, .(cell, cluster, major_celltype, epithelial, immune, stromal)],
       file.path(OUT, sprintf("%s_res%s_percell_module_scores.csv", SAMPLE, RES)))

print(cmp[, .(cluster, n_cells, old_call, ams_call, agrees,
              old_margin, ams_margin, old_ambiguous, ams_ambiguous)])
disagree <- cmp[agrees == FALSE]
log(sprintf("clusters disagreeing: %d of %d (%.1f%% of cells)", nrow(disagree), nrow(cmp),
            100 * sum(disagree$n_cells) / sum(cmp$n_cells)))
log(sprintf("ambiguous: frozen %d, AddModuleScore %d",
            sum(cmp$old_ambiguous), sum(cmp$ams_ambiguous)))
log("done -> ", OUT)
