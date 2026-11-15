#!/usr/bin/env Rscript
# res0.3 clusters on the tissue, for inspection.
#
# Two figures, because one map cannot do both jobs at 715,479 cells:
#   *_clusters_spatial.{pdf,png}   all twelve at once, high-contrast palette
#   *_clusters_facets.png          one panel per cluster, that cluster against grey
#
# Cluster sizes run from 210,750 cells (29%) down to 641 (0.09%), so on the combined
# map the small clusters are drawn last and a little larger -- otherwise cluster 11
# is 641 points under half a million and simply is not there. The facet grid is the
# one to trust for the small clusters.
#
# Palette: twelve hues picked to stay apart on white at a 0.1 pt dot, with the
# similar pairs (the two blues, the two reds) assigned to clusters that do not touch
# on the tissue. The atera palette the supplement uses only separates cleanly up to
# about eight.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05f_cluster_spatial.R

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

d <- fread(file.path(CLUST, sprintf("%s_annotated_cells_res%s.csv", SAMPLE, RES)))
setnames(d, "V1", "cell")
ann <- fread(file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES)))
ann[, cluster := as.character(cluster)]

d[, cluster := as.character(cluster_frozen)]
lv <- as.character(sort(as.integer(ann$cluster)))
d[, cluster := factor(cluster, levels = lv)]
d[, y := -y]                       # image y grows downward, as in Fig 5's panels
log("cells: ", nrow(d), "  clusters: ", nlevels(d$cluster))

PAL <- c("0"  = "#E6194B", "1"  = "#3CB44B", "2"  = "#4363D8", "3"  = "#F58231",
         "4"  = "#911EB4", "5"  = "#008080", "6"  = "#9A6324", "7"  = "#F032E6",
         "8"  = "#000075", "9"  = "#808000", "10" = "#00CED1", "11" = "#800000")
PAL <- PAL[lv]

n_of <- setNames(ann$n_cells[match(lv, ann$cluster)], lv)
ct   <- setNames(ann$major_celltype[match(lv, ann$cluster)], lv)
amb  <- setNames(ann$ambiguous[match(lv, ann$cluster)], lv)
lab  <- sprintf("%s  %s%s  (%s)", lv, ct[lv], ifelse(amb[lv], "*", ""),
                formatC(n_of[lv], big.mark = ",", format = "d"))
names(lab) <- lv

th <- theme_void(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 8, face = "bold", hjust = .5, margin = margin(b = 2)),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.size = unit(0.34, "cm"))

## ---- 1. all clusters on one map ---------------------------------------------
# smallest last so it lands on top, and sized by rank so 641 cells stay visible
sz <- setNames(seq(0.055, 0.5, length.out = length(lv))[rank(n_of[lv])], lv)
d[, psize := sz[as.character(cluster)]]
setorder(d, -psize)                 # large points (small clusters) drawn last
p1 <- ggplot(d, aes(x, y, colour = cluster)) +
  geom_point(aes(size = psize), stroke = 0, show.legend = c(colour = TRUE, size = FALSE)) +
  scale_colour_manual(values = PAL, labels = lab, name = "Cluster (res 0.3)") +
  scale_size_identity() + coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2.4), ncol = 1)) +
  labs(title = sprintf("%s res %s -- all %d clusters", SAMPLE, RES, length(lv))) + th

asp <- diff(range(d$y)) / diff(range(d$x))
W1 <- 9; H1 <- (W1 - 2.1) * asp + 0.5
for (ext in c("pdf", "png")) {
  f <- file.path(OUT, sprintf("%s_res%s_clusters_spatial.%s", SAMPLE, RES, ext))
  if (ext == "png" && requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(f, width = W1, height = H1, units = "in", res = 400, background = "white")
  } else if (ext == "png") {
    png(f, width = W1, height = H1, units = "in", res = 400, bg = "white", type = "cairo")
  } else {
    grDevices::cairo_pdf(f, width = W1, height = H1)
  }
  print(p1); invisible(dev.off())
}
log(sprintf("wrote clusters_spatial  (%.1f x %.1f in)", W1, H1))

## ---- 2. one panel per cluster ------------------------------------------------
bg <- d[, .(x, y)]
facet1 <- function(g) {
  fg <- d[cluster == g]
  ggplot() +
    geom_point(data = bg, aes(x, y), colour = "grey90", size = 0.04, stroke = 0) +
    geom_point(data = fg, aes(x, y), colour = PAL[[g]], size = 0.07, stroke = 0) +
    coord_fixed() + labs(title = lab[[g]]) + th +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = .5, colour = PAL[[g]]))
}
p2 <- wrap_plots(lapply(lv, facet1), ncol = 4)
W2 <- 9; H2 <- 3 * ((W2 / 4) * asp + 0.22)
f2 <- file.path(OUT, sprintf("%s_res%s_clusters_facets.png", SAMPLE, RES))
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(f2, width = W2, height = H2, units = "in", res = 400, background = "white")
} else {
  png(f2, width = W2, height = H2, units = "in", res = 400, bg = "white", type = "cairo")
}
print(p2); invisible(dev.off())
log(sprintf("wrote clusters_facets  (%.1f x %.1f in)", W2, H2))
log("done -> ", OUT)
