#!/usr/bin/env Rscript
# All three compartments scored against each other, per cell, with AddModuleScore.
#
# novaseq_05k asked the one-sided question -- do epithelial markers separate epi from
# the rest -- and the answer was no, they separate immune from the rest. This asks the
# symmetric version: score every cell for all three modules and see whether its own
# compartment's module wins.
#
# AddModuleScore rather than the pipeline's z-score because it is absolute: a cell's
# score is its module mean minus the mean of ~100 expression-matched control genes, so
# a cell that expresses none of the module lands near 0 instead of being ranked against
# the other groups. The pipeline's recipe z-scores across clusters, which forces one
# cluster to win every module however little of it anyone expresses.
#
# Panels are the three frozen compartments only -- the per-cluster version is in 05k.
#
# The panels are the classic ones and deliberately larger than the five-gene sets the
# pipeline uses: AddModuleScore's control correction is noisy with five genes, and
# MUC1 (one of the pipeline's five) is expressed in 0.5% of the epithelial cells here.
# HLA-DRA and CD74 are left out of the immune panel on purpose -- breast epithelium
# expresses both.
#
# Run (~1 min; reads the capped matrix, not the 7 GB object):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05l_three_module_addmodulescore.R

suppressWarnings(suppressMessages({
  library(data.table); library(Matrix); library(ggplot2); library(patchwork)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")

MARKERS <- list(
  epithelial = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7",
                 "CLDN3", "CLDN4", "CLDN7", "ELF3", "KRT5", "KRT14", "KRT17"),
  immune     = c("PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "TRBC2", "IL7R",
                 "CD68", "CD14", "LYZ", "TYROBP", "AIF1", "C1QA", "C1QB", "FCER1G",
                 "MS4A1", "CD79A", "JCHAIN", "MZB1", "TPSAB1", "CPA3", "MS4A2"),
  stromal    = c("COL1A1", "COL1A2", "COL3A1", "COL6A2", "DCN", "LUM", "FBLN1",
                 "POSTN", "FAP", "PDGFRB", "ACTA2", "MCAM", "RGS5",
                 "VWF", "PECAM1", "CDH5", "CLDN5"))
MAP    <- c(epithelial = "epi", immune = "imm", stromal = "stroma")
CMP_LV <- c("epi", "imm", "stroma")
PAL    <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

b   <- readRDS(file.path(OUT, sprintf("%s_res%s_markers_input.rds", SAMPLE, RES)))
ann <- fread(file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES)))
ann[, cluster := as.character(cluster)]
froz <- setNames(ann$major_celltype, ann$cluster)
X <- b$data
cmp <- factor(froz[as.character(b$idents)], levels = CMP_LV)
log("cells: ", ncol(X), " (2,000 per cluster)  genes: ", nrow(X))
for (k in names(MARKERS)) {
  miss <- setdiff(MARKERS[[k]], rownames(X))
  if (length(miss)) log("  ", k, ": not in the matrix, dropped -- ", paste(miss, collapse = " "))
  MARKERS[[k]] <- intersect(MARKERS[[k]], rownames(X))
}

# Seurat 5.1.0's AddModuleScore body, reproduced in novaseq_05b because SeuratObject
# 5.4.0 made GetAssayData(slot=) defunct; here it takes the matrix directly.
add_module_score <- function(mat, features, nbin = 24, ctrl = 100, seed = 1) {
  set.seed(seed)
  data.avg <- Matrix::rowMeans(mat)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(data.avg + rnorm(length(data.avg)) / 1e30,
                                  n = nbin, labels = FALSE, right = FALSE)
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector("list", length(features))
  for (i in seq_along(features))
    for (g in features[[i]])
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(sample(data.cut[which(data.cut == data.cut[g])],
                                      size = ctrl, replace = FALSE)))
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- t(sapply(ctrl.use, function(g) Matrix::colMeans(mat[g, ])))
  feat.scores <- t(sapply(features, function(g) Matrix::colMeans(mat[g, , drop = FALSE])))
  out <- as.data.table(t(feat.scores - ctrl.scores))
  setnames(out, names(features)); out
}

log("AddModuleScore (nbin 24, ctrl 100, seed 1)")
s <- add_module_score(X, MARKERS)
s[, `:=`(cell = colnames(X), frozen = cmp)]
s[, ams := factor(MAP[names(MARKERS)[max.col(as.matrix(.SD))]], levels = CMP_LV),
  .SDcols = names(MARKERS)]
fwrite(s, file.path(OUT, "strataMap_percell_three_module_ams.csv"))

cat("\n=== median score per module, by frozen compartment ===\n")
m <- s[, lapply(.SD, function(z) round(median(z), 3)), by = frozen, .SDcols = names(MARKERS)]
print(m)
cat("\n=== % of cells with a positive score ===\n")
print(s[, lapply(.SD, function(z) round(100 * mean(z > 0), 1)), by = frozen,
        .SDcols = names(MARKERS)])
cat("\n=== argmax of the three scores (rows frozen, cols AddModuleScore) ===\n")
tb <- table(s$frozen, s$ams); print(addmargins(tb))
cat(sprintf("\nagreement %.1f%% over %d cells (2,000 per cluster, so compartments are\n",
            100 * sum(diag(tb)) / sum(tb), nrow(s)))
cat("equally weighted here and not in proportion to their size)\n")

th <- theme_bw(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        legend.key.size = unit(0.26, "cm"), legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7))

## ---- a: the three scores, by frozen compartment ------------------------------
lng <- melt(s, id.vars = c("cell", "frozen"), measure.vars = names(MARKERS),
            variable.name = "module", value.name = "score")
lng[, module := factor(MAP[as.character(module)], levels = CMP_LV)]
lng[, mlab := factor(sprintf("%s module\n(%d genes)", CMP_LV,
                             lengths(MARKERS)[match(CMP_LV, MAP)])[match(module, CMP_LV)],
                     levels = sprintf("%s module\n(%d genes)", CMP_LV,
                                      lengths(MARKERS)[match(CMP_LV, MAP)]))]
p_a <- ggplot(lng, aes(frozen, score, fill = frozen)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey45") +
  geom_violin(scale = "width", linewidth = .2, colour = "grey30") +
  geom_boxplot(width = .16, outlier.shape = NA, linewidth = .22, fill = "white") +
  facet_wrap(~ mlab, nrow = 1) +
  scale_fill_manual(values = PAL, name = "frozen call") +
  labs(title = "a  All three modules, scored per cell",
       subtitle = "AddModuleScore: module mean minus ~100 expression-matched control genes, so 0 means no enrichment",
       x = "res0.3 frozen compartment", y = "module score") + th
log("violins done")

## ---- b: where the argmax lands -----------------------------------------------
cf <- as.data.table(tb); setnames(cf, c("frozen", "ams", "n"))
cf[, pct := 100 * n / sum(n), by = frozen]
cf[, `:=`(frozen = factor(frozen, levels = rev(CMP_LV)),
          ams = factor(ams, levels = CMP_LV))]
p_b <- ggplot(cf, aes(ams, frozen, fill = pct)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.1f%%\n%s", pct, formatC(n, big.mark = ",", format = "d"))),
            size = 2.3, colour = ifelse(cf$pct > 55, "white", "grey20")) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B"),
                       limits = c(0, 100), name = "% of row") +
  coord_fixed() +
  labs(title = "b  Which module wins",
       subtitle = "row = frozen call, column = argmax of the three scores",
       x = "AddModuleScore argmax", y = "frozen compartment") + th

## ---- c: the margin ------------------------------------------------------------
M <- as.matrix(s[, names(MARKERS), with = FALSE])
srt <- t(apply(M, 1L, sort, decreasing = TRUE))
s[, margin := srt[, 1] - srt[, 2]]
p_c <- ggplot(s, aes(margin, colour = frozen)) +
  stat_ecdf(linewidth = .5) +
  scale_colour_manual(values = PAL, name = "frozen call") +
  labs(title = "c  How clear is the winner",
       subtitle = "ECDF of (best score - second best); left means the top two modules are tied",
       x = "margin", y = "fraction of cells") + th
cat("\n=== margin between the best and second-best module ===\n")
print(s[, .(n = .N, median_margin = round(median(margin), 3),
            pct_margin_lt_0.1 = round(100 * mean(margin < 0.1), 1)), by = frozen])

fig <- p_a / (p_b | p_c) + plot_layout(heights = c(1, 1))
W <- 8.3; H <- 6.4
for (ext in c("pdf", "png")) {
  f <- file.path(OUT, sprintf("%s_res%s_three_module_ams.%s", SAMPLE, RES, ext))
  if (ext == "png" && requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(f, width = W, height = H, units = "in", res = 400, background = "white")
  } else if (ext == "png") {
    png(f, width = W, height = H, units = "in", res = 400, bg = "white", type = "cairo")
  } else {
    grDevices::cairo_pdf(f, width = W, height = H)
  }
  print(fig); invisible(dev.off())
  log("wrote ", basename(f))
}
log("done -> ", OUT)
