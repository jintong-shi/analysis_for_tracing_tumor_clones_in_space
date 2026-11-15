#!/usr/bin/env Rscript
# Is the Ig signal in this sample made by the cells that carry it, or is it ambient?
#
# Six panels, each a different way of asking the same question. Reads the per-cell Ig
# counts novaseq_05i wrote and the res0.3 / res1 calls.
#
#  a  Ig fraction on the tissue              a field or discrete cells?
#  b  JCHAIN fraction on the tissue           the control: plasma machinery, not secreted
#                                             in bulk, so it contaminates far less
#  c  per res1 cluster, mean Ig vs mean       a real plasma cluster is high in both; a
#     JCHAIN                                  contaminated one is high in Ig alone
#  d  Ig fraction by res1 major type          endothelium/adipocyte CANNOT make Ig, so
#                                             whatever they carry is the ambient floor
#  e  Ig in non-plasma cells vs distance      ambient diffuses, so it must decay away
#     from the nearest plasma-rich bin        from the source; cell-intrinsic Ig must not
#  f  per 100 um bin: local plasma density    the same argument without a distance
#     vs Ig in that bin's NON-plasma cells    threshold
#
# MZB1 was the first choice for b and c and had to be dropped: only 34,808 counts in the
# whole sample and 3.1% of cells positive, too sparse to say anything per cell. JCHAIN
# has 480,157 counts and 13% positive, and separates plasma from tumour 54-fold
# (4.66 vs 0.09 counts per cell) against Ig's 14-fold.
#
# Light-chain restriction (kappa vs lambda per cell) was the fourth test and is NOT
# shown, because it does not discriminate here: among cells with >=50 light-chain UMI,
# 3.1% of the plasma cells are lambda-dominant against 2.1% of the disputed ones, and
# the plasma median kappa share is 0.90 against a whole-tissue 0.83. Either this
# infiltrate is near-clonally kappa, or even the plasma cluster's light chains are
# ambient-dominated. The data cannot separate those, so the panel would mislead.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05j_ig_ambient_figs.R

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
PT     <- 0.07
BIN    <- 100        # um; the StrataMap grid is in um, so a bin is 100 x 100 um

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

ig <- fread(file.path(OUT, sprintf("%s_percell_ig_counts.csv", SAMPLE)))
b  <- fread(file.path(OUT, sprintf("%s_res1_cells_with_major.csv", SAMPLE)))

d <- merge(b[, .(cell, x, y, c1 = cluster, maj1 = major)], ig, by = "cell")
stopifnot(nrow(d) == nrow(b))

d[, ig_frac := 100 * (IGH_C + IGK + IGL) / nCount]
d[, jc_frac := 100 * JCHAIN / nCount]
log("cells: ", nrow(d), sprintf("  median Ig %.2f%%  max %.1f%%",
                                median(d$ig_frac), max(d$ig_frac)))

# Clusters are named by number, not by the label in the assignment table, because this
# figure is what moved cluster 4 OUT of Plasma cell: keyed on the label it would quietly
# change meaning the moment the table was corrected.
#   13  the plasma cluster both resolutions agree on
#    4  44,877 cells first called plasma off 10/10 Ig top genes, now tumour
CL_PLASMA <- 13L; CL_DISPUTED <- 4L
d[, grp := fifelse(c1 == CL_PLASMA, "plasma (cluster 13)",
           fifelse(c1 == CL_DISPUTED, "disputed (cluster 4)",
           fifelse(maj1 %in% c("Endothelial", "Adipocyte", "Basal / myoepithelial",
                               "Pericyte / SMC"), "cannot make Ig", "other")))]
GRP_LV  <- c("plasma (cluster 13)", "disputed (cluster 4)", "cannot make Ig", "other")
GRP_PAL <- c("#B2182B", "#E69F00", "#2166AC", "grey70"); names(GRP_PAL) <- GRP_LV
d[, grp := factor(grp, levels = GRP_LV)]
print(d[, .(n = .N, median_ig = round(median(ig_frac), 2),
            mean_jchain = round(mean(JCHAIN), 3)), by = grp][order(-n)])

th <- theme_bw(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        legend.key.size = unit(0.26, "cm"), legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7))
thv <- theme_void(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 8, face = "bold", hjust = .5, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, hjust = .5, colour = "grey30",
                                     margin = margin(b = 2)),
        legend.position = "bottom", legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5), legend.key.size = unit(0.26, "cm"))

## ---- a, b: the two fields ----------------------------------------------------
sp <- function(v, ttl, sub, hi) {
  dd <- d[order(get(v))]                      # high values last so they are visible
  ggplot(dd, aes(x, y, colour = pmin(get(v), hi))) +
    geom_point(size = PT, stroke = 0) + coord_fixed() +
    scale_colour_gradientn(colours = c("grey92", "#FEE090", "#F46D43", "#A50026"),
                           limits = c(0, hi), name = "% of UMI") +
    guides(colour = guide_colorbar(barwidth = unit(2, "cm"), barheight = unit(0.22, "cm"),
                                   title.position = "top")) +
    labs(title = ttl, subtitle = sub) + thv
}
p_a <- sp("ig_frac", "a  Immunoglobulin",
          "IGH-C + IGK + IGL, capped at 10%", 10)
p_b <- sp("jc_frac", "b  JCHAIN (plasma-cell machinery)",
          "not secreted in bulk, capped at 1%", 1)

## ---- c: one point per res1 cluster -------------------------------------------
# per cell JCHAIN is too sparse to scatter (13% of cells positive), but a cluster mean
# over 2,000+ cells is stable, and it is the cluster that got the label anyway
cl <- d[, .(n = .N, ig = mean(ig_frac), jc = mean(jc_frac),
            maj = maj1[1]), by = c1]
cl[, lab := fifelse(maj == "Plasma cell" | ig > 1, as.character(c1), NA_character_)]
p_c <- ggplot(cl, aes(ig, jc, size = n, colour = maj == "Plasma cell")) +
  geom_point(alpha = .75, stroke = 0) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.1, colour = "grey20",
                           min.segment.length = 0, segment.size = .2, show.legend = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_size_continuous(range = c(1, 4.5), name = "cells") +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "called plasma", `FALSE` = "other"), name = NULL) +
  labs(title = "c  Ig product vs plasma machinery, per res1 cluster",
       subtitle = "a true plasma cluster is high in both; ambient lifts Ig only",
       x = "mean Ig % of UMI", y = "mean JCHAIN % of UMI") + th

## ---- d: the ambient floor ----------------------------------------------------
ord <- d[, .(m = median(ig_frac)), by = maj1][order(-m)]$maj1
d[, maj1f := factor(maj1, levels = rev(ord))]
cant <- d[grp == "cannot make Ig", median(ig_frac)]
p_d <- ggplot(d, aes(maj1f, ig_frac + 0.01)) +
  geom_hline(yintercept = cant + 0.01, linetype = 2, colour = "#2166AC", linewidth = .35) +
  geom_boxplot(outlier.shape = NA, linewidth = .25, fill = "grey93", width = .62) +
  scale_y_log10() + coord_flip(ylim = c(0.03, 100)) +
  labs(title = "d  Ig by res1 major type",
       subtitle = sprintf("dashed = median of the lineages that cannot make Ig (%.2f%%)", cant),
       x = NULL, y = "Ig % of UMI") + th

## ---- e: ambient decays with distance from the source -------------------------
d[, `:=`(bx = floor(x / BIN), by = floor(y / BIN))]
gb <- d[, .(n = .N, plasma_pct = 100 * mean(maj1 == "Plasma cell"),
            ig_nonplasma = median(ig_frac[maj1 != "Plasma cell"])),
        by = .(bx, by)][n >= 30 & !is.na(ig_nonplasma)]
src <- as.matrix(gb[plasma_pct >= 50, .(bx, by)])          # the plasma-rich bins
stopifnot(nrow(src) > 0)
gb[, dist_um := BIN * apply(as.matrix(.SD), 1L, function(p)
     sqrt(min((src[, 1] - p[1])^2 + (src[, 2] - p[2])^2))), .SDcols = c("bx", "by")]
dec <- gb[plasma_pct < 50][, db := pmin(floor(dist_um / 200), 10L)][
  , .(n = .N, m = median(ig_nonplasma),
      lo = quantile(ig_nonplasma, .25), hi = quantile(ig_nonplasma, .75)), by = db][order(db)]
dec[, mid := db * 200 + 100]
flo <- d[grp == "cannot make Ig", median(ig_frac)]
p_e <- ggplot(dec, aes(mid, m)) +
  geom_hline(yintercept = flo, linetype = 2, colour = "#2166AC", linewidth = .35) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#F46D43", alpha = .2) +
  geom_line(colour = "#A50026", linewidth = .5) +
  geom_point(colour = "#A50026", size = 1) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 2000, 400),
                     labels = c(seq(0, 1800, 400), ">2000")[c(1:5, 6)]) +
  labs(title = "e  Ambient decays away from the plasma foci",
       subtitle = "non-plasma cells only, IQR shaded; dashed = floor",
       x = "distance to the nearest plasma-rich bin (um)", y = "Ig % of UMI") + th
log("decay: ", paste(sprintf("%dum %.2f%%", dec$mid, dec$m), collapse = "  "))

## ---- f: the same without a distance threshold --------------------------------
ct <- cor(gb$plasma_pct, gb$ig_nonplasma, method = "spearman")
p_f <- ggplot(gb, aes(plasma_pct, ig_nonplasma + 0.01)) +
  geom_point(size = .3, stroke = 0, colour = "#4C72B0", alpha = .4) +
  geom_smooth(method = "loess", se = FALSE, colour = "#B2182B", linewidth = .6, span = .6) +
  scale_y_log10() + scale_x_sqrt(breaks = c(0, 1, 5, 15, 40, 100)) +
  labs(title = sprintf("f  Ambient tracks local plasma density (rho %.2f)", ct),
       subtitle = sprintf("%s bins of %d um with >=30 cells; y is that bin's non-plasma cells",
                          format(nrow(gb), big.mark = ","), BIN),
       x = "% plasma cells in the bin (sqrt scale)", y = "Ig % of UMI, non-plasma cells") + th
log(sprintf("bins: %d  spearman rho = %.3f", nrow(gb), ct))

fig <- (p_a | p_b | p_c) / (p_d | p_e | p_f)
W <- 11; H <- 8.2
for (ext in c("pdf", "png")) {
  f <- file.path(OUT, sprintf("%s_ig_ambient.%s", SAMPLE, ext))
  if (ext == "png" && requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(f, width = W, height = H, units = "in", res = 400, background = "white")
  } else if (ext == "png") {
    png(f, width = W, height = H, units = "in", res = 400, bg = "white", type = "cairo")
  } else {
    grDevices::cairo_pdf(f, width = W, height = H)
  }
  suppressWarnings(print(fig)); invisible(dev.off())
  log("wrote ", basename(f))
}
log("done -> ", OUT)
