#!/usr/bin/env Rscript
# Where the res0.3 and the res1 top-gene calls disagree, on the tissue.
#
# Both calls come from FindAllMarkers top-10 genes read by hand -- res0.3 from
# novaseq_05d (12 clusters), res1 from the same script at res 1 (31 clusters). The
# assignment tables are the two *_major_assignment_with_compartment.csv files.
#
# Row 1  the two epi/rest calls side by side, then the cells that changed
# Row 2  the five largest reclassification flows, one panel each
#
# Every point is PT across every panel: a dot's area means one cell and nothing
# else. Only draw order is tuned -- ggplot paints in row order, so the changed
# cells go last and are not buried under the 611,045 that agree.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05h_res03_vs_res1_spatial.R

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
PT     <- 0.07
N_FLOW <- 5L

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## ---- join the two calls cell by cell ----------------------------------------
a <- fread(file.path(CLUST, sprintf("%s_annotated_cells_res0.3.csv", SAMPLE)))
setnames(a, "V1", "cell")
a[, c03 := as.character(cluster_frozen)]
m0 <- fread(file.path(OUT, sprintf("%s_res0.3_major_assignment_with_compartment.csv", SAMPLE)))
m0[, c03 := as.character(cluster)]
a[m0, on = "c03", `:=`(maj0 = i.major, epi0 = i.epi)]

b <- fread(file.path(OUT, sprintf("%s_res1_cells_with_major.csv", SAMPLE)))
d <- merge(a[, .(cell, c03, maj0, epi0)],
           b[, .(cell, x, y, maj1 = major, epi1 = epi)], by = "cell")   # x/y already flipped
stopifnot(nrow(d) == nrow(b))
log("cells: ", nrow(d), "  disagreeing: ", sum(d$epi0 != d$epi1))

EPI_LV  <- c("Epithelium", "Rest")
EPI_PAL <- c(Epithelium = "#C44E52", Rest = "#4C72B0")
d[, epi0 := factor(epi0, levels = EPI_LV)]
d[, epi1 := factor(epi1, levels = EPI_LV)]
d[, chg := fifelse(epi0 == epi1, "unchanged",
                   fifelse(epi0 == "Epithelium", "Epi (res0.3) -> Rest (res1)",
                                                 "Rest (res0.3) -> Epi (res1)"))]
CHG_LV  <- c("Epi (res0.3) -> Rest (res1)", "Rest (res0.3) -> Epi (res1)", "unchanged")
CHG_PAL <- c("#E69F00", "#009E73", "grey88"); names(CHG_PAL) <- CHG_LV
d[, chg := factor(chg, levels = CHG_LV)]

asp <- diff(range(d$y)) / diff(range(d$x))
th <- theme_void(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 8, face = "bold", hjust = .5, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, hjust = .5, colour = "grey30",
                                     margin = margin(b = 2)),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5), legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 2, 2, 2))

n_fmt <- function(n) formatC(n, big.mark = ",", format = "d")

# rare level last so it lands on top; sizes are never touched
spat <- function(col, pal, ttl, sub) {
  o  <- names(sort(table(d[[col]]), decreasing = TRUE))
  dd <- d[order(match(as.character(get(col)), o), decreasing = TRUE)]
  ggplot(dd, aes(x, y, colour = .data[[col]])) +
    geom_point(size = PT, stroke = 0) + coord_fixed() +
    scale_colour_manual(values = pal, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 1.8), nrow = 2)) +
    labs(title = ttl, subtitle = sub) + th
}

## ---- row 1 -------------------------------------------------------------------
p1 <- spat("epi0", EPI_PAL, "res 0.3 call (12 clusters)",
           sprintf("Epithelium %s  |  Rest %s", n_fmt(sum(d$epi0 == "Epithelium")),
                   n_fmt(sum(d$epi0 == "Rest"))))
p2 <- spat("epi1", EPI_PAL, "res 1 call (31 clusters)",
           sprintf("Epithelium %s  |  Rest %s", n_fmt(sum(d$epi1 == "Epithelium")),
                   n_fmt(sum(d$epi1 == "Rest"))))
p3 <- spat("chg", CHG_PAL, "cells whose call changed",
           sprintf("%s cells differ (%.2f%% agreement)",
                   n_fmt(sum(d$epi0 != d$epi1)), 100 * mean(d$epi0 == d$epi1)))

## ---- row 2: the largest major-type flows ------------------------------------
fl <- d[maj0 != maj1, .N, by = .(maj0, maj1)][order(-N)][seq_len(N_FLOW)]
bg <- d[, .(x, y)]
flow1 <- function(i) {
  fg <- d[maj0 == fl$maj0[i] & maj1 == fl$maj1[i]]
  ggplot() +
    geom_point(data = bg, aes(x, y), colour = "grey91", size = PT, stroke = 0) +
    geom_point(data = fg, aes(x, y), colour = "#D55E00", size = PT, stroke = 0) +
    coord_fixed() +
    labs(title = sprintf("%s  ->  %s", fl$maj0[i], fl$maj1[i]),
         subtitle = sprintf("%s cells", n_fmt(fl$N[i]))) +
    th + theme(legend.position = "none",
               plot.title = element_text(size = 7, face = "bold", hjust = .5,
                                         margin = margin(b = 1)))
}
log("top flows: ", paste(sprintf("%s->%s (%d)", fl$maj0, fl$maj1, fl$N), collapse = "; "))

fig <- wrap_plots(c(list(p1, p2, p3), lapply(seq_len(N_FLOW), flow1)),
                  ncol = 4, heights = c(1, 1))

W <- 11; H <- 2 * ((W / 4) * asp + 0.85)
for (ext in c("pdf", "png")) {
  f <- file.path(OUT, sprintf("%s_res0.3_vs_res1_discordance_spatial.%s", SAMPLE, ext))
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
log(sprintf("done (%.1f x %.1f in) -> %s", W, H, OUT))
