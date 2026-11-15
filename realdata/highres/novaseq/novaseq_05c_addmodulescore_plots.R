#!/usr/bin/env Rscript
# Check sheet for novaseq_05b_addmodulescore_check.R -- the AddModuleScore result
# put on the UMAP and on the tissue so the three flipped clusters can be eyeballed.
#
# Rows:
#   1  epithelial / immune / stromal per-cell score on the UMAP
#   2  the same three on the tissue
#   3  frozen compartment call | AddModuleScore call | the cells that changed, UMAP
#   4  the same three on the tissue
#
# The scores are per cell; the CALL is per cluster (mean score, then argmax), which
# is why rows 3-4 are flat within a cluster while rows 1-2 are not.
#
# PNG only: 715,479 cells x 12 panels as vector runs to a few hundred MB and this is
# a check sheet, not a figure.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05c_addmodulescore_plots.R

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")

MODS <- c(epithelial = "Epithelial", immune = "Immune", stromal = "Stromal")
PAL  <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
LEV  <- c("epi", "imm", "stroma")
BASE <- 7

sc  <- fread(file.path(OUT, sprintf("%s_res%s_percell_module_scores.csv", SAMPLE, RES)))
xy  <- fread(file.path(CLUST, sprintf("%s_annotated_cells_res%s.csv", SAMPLE, RES)))
setnames(xy, "V1", "cell")
cmp <- fread(file.path(OUT, sprintf("%s_res%s_addmodulescore_vs_frozen.csv", SAMPLE, RES)))
cmp[, cluster := as.character(cluster)]

d <- merge(sc[, .(cell, cluster = as.character(cluster), epithelial, immune, stromal)],
           xy[, .(cell, x, y, umap1, umap2)], by = "cell")
d[cmp, `:=`(frozen = i.old_call, ams = i.ams_call, flipped = !i.agrees), on = "cluster"]
# image y grows downward, as in Fig 5's StrataMap panels
d[, y := -y]
d[, frozen := factor(frozen, levels = LEV)]
d[, ams    := factor(ams,    levels = LEV)]
d[, change := factor(fifelse(flipped, paste0(frozen, " -> ", ams), "unchanged"))]
stopifnot(nrow(d) == nrow(sc))
message(sprintf("cells: %d  flipped: %d (%.1f%%)", nrow(d), sum(d$flipped),
                100 * mean(d$flipped)))

chg_lv  <- c(setdiff(levels(d$change), "unchanged"), "unchanged")
d[, change := factor(change, levels = chg_lv)]
chg_pal <- setNames(c("#E69F00", "#0072B2", "#CC79A7")[seq_len(length(chg_lv) - 1)],
                    head(chg_lv, -1))
chg_pal <- c(chg_pal, unchanged = "grey88")

th <- theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = .5, margin = margin(b = 1)),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
        legend.key.size = unit(0.26, "cm"), legend.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(2, 2, 2, 2))

# diverging scale centred on 0: AddModuleScore is a difference from matched controls,
# so the sign is what matters and a symmetric limit keeps the three panels comparable
lim <- max(abs(quantile(unlist(d[, names(MODS), with = FALSE]), c(.005, .995))))
feat <- function(m, xv, yv, fixed) {
  dd <- d[order(abs(get(m)))]
  ggplot(dd, aes(.data[[xv]], .data[[yv]], colour = .data[[m]])) +
    geom_point(size = .05, stroke = 0) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0,
                           limits = c(-lim, lim), oob = scales::squish, name = "Score") +
    guides(colour = guide_colorbar(barwidth = unit(1.5, "cm"), barheight = unit(0.22, "cm"),
                                   title.position = "top")) +
    labs(title = MODS[[m]]) + (if (fixed) coord_fixed() else NULL) + th
}
disc <- function(v, pal, nm, ttl, xv, yv, fixed) {
  dd <- if (v == "change") rbind(d[change == "unchanged"], d[change != "unchanged"]) else d
  ggplot(dd, aes(.data[[xv]], .data[[yv]], colour = .data[[v]])) +
    geom_point(size = .05, stroke = 0) +
    scale_colour_manual(values = pal, name = nm, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 1.6), nrow = 2, title.position = "top")) +
    labs(title = ttl) + (if (fixed) coord_fixed() else NULL) + th
}

mods <- names(MODS)
r1 <- lapply(mods, feat, xv = "umap1", yv = "umap2", fixed = FALSE)
r2 <- lapply(mods, feat, xv = "x",     yv = "y",     fixed = TRUE)
r3 <- list(disc("frozen", PAL, "Compartment", "Frozen call (marker z-score)", "umap1", "umap2", FALSE),
           disc("ams",    PAL, "Compartment", "AddModuleScore call",          "umap1", "umap2", FALSE),
           disc("change", chg_pal, "Change",  "Cells whose call changed",     "umap1", "umap2", FALSE))
r4 <- list(disc("frozen", PAL, "Compartment", "Frozen call (marker z-score)", "x", "y", TRUE),
           disc("ams",    PAL, "Compartment", "AddModuleScore call",          "x", "y", TRUE),
           disc("change", chg_pal, "Change",  "Cells whose call changed",     "x", "y", TRUE))

asp <- diff(range(d$y)) / diff(range(d$x))
h_um <- 2.3; h_sp <- (7.09 / 3) * asp + 0.7
W <- 7.09; H <- 2 * h_um + 2 * h_sp + 0.2
fig <- wrap_plots(c(r1, r2, r3, r4), ncol = 3,
                  heights = c(h_um, h_sp, h_um, h_sp)[c(1, 2, 3, 4)])
f <- file.path(OUT, sprintf("%s_res%s_addmodulescore_check.png", SAMPLE, RES))
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(f, width = W, height = H, units = "in", res = 300, background = "white")
} else {
  png(f, width = W, height = H, units = "in", res = 300, bg = "white", type = "cairo")
}
print(fig); invisible(dev.off())
message(sprintf("wrote %s  (%.1f x %.1f in)", f, W, H))
