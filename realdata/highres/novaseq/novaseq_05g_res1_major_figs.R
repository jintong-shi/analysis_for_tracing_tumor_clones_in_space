#!/usr/bin/env Rscript
# res1 clusters and the major cell types read off their markers.
#
# Inputs, all already on disk:
#   <CLUST>/IDC_Grade3_cell_metadata.csv            x, y, umap_1, umap_2, spatial_snn_res.1
#   <OUT>/IDC_Grade3_res1_top10_markers.csv         novaseq_05d step 'test'
#   <OUT>/IDC_Grade3_res1_markers_input.rds         novaseq_05d step 'extract' (2,000 cells/cluster)
#   <OUT>/IDC_Grade3_res1_major_assignment.csv      cluster -> major, written by hand from
#                                                   the top markers; this script never guesses
#
# Outputs, all in <OUT>:
#   *_res1_clusters_spatial.png     all res1 clusters, high-contrast
#   *_res1_clusters_facets.png      one panel per cluster
#   *_res1_top_marker_dotplot.pdf   top markers x cluster
#   *_res1_major_spatial_umap.png   predicted major types, tissue + UMAP
#   *_res1_major_classic_dotplot.pdf classic lineage markers x predicted major type
#
# Expression for both dot plots comes from the capped matrix, so the 7 GB object is
# never reopened. Dot size is the percent of cells expressing and dot colour is the
# mean of the expressing-and-not cells, z-scored across groups -- the usual pair, and
# the two say different things, so neither is dropped.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05g_res1_major_figs.R

suppressWarnings(suppressMessages({
  library(data.table); library(Matrix); library(ggplot2); library(patchwork)
}))

PT <- 0.07          # one point size for every spatial and UMAP panel, so that what
                    # the eye reads as area is abundance and nothing else
PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "1"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
TOP_N  <- 5L           # per cluster in the dot plot; the table keeps 10

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
st <- function(f, ...) file.path(OUT, sprintf(paste0("%s_res%s_", f), SAMPLE, RES, ...))

# canonical markers per lineage, for the major-type dot plot. Deliberately NOT the
# genes the clusters were called on, so the panel is a check and not a restatement.
CLASSIC <- list(
  `Tumour / luminal` = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "AGR2", "CPB1"),
  `Basal / myoepi`   = c("KRT5", "KRT14", "KRT17", "MIR205HG", "TP63", "ACTA2"),
  `Plasma`           = c("MZB1", "JCHAIN", "IGKC", "IGHG1", "DERL3"),
  `B cell`           = c("MS4A1", "CD79A", "BANK1"),
  `T / NK`           = c("CD3D", "CD3E", "CD2", "IL7R", "NKG7"),
  `Myeloid`          = c("CD68", "CD14", "LYZ", "C1QA", "ITGAX"),
  `Mast`             = c("TPSAB1", "CPA3", "MS4A2"),
  `Fibroblast / CAF` = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1", "POSTN"),
  `Endothelial`      = c("PECAM1", "VWF", "CLDN5", "CDH5"),
  `Pericyte / SMC`   = c("RGS5", "MYH11", "NOTCH3", "TAGLN"),
  `Adipocyte`        = c("ADIPOQ", "FABP4"),
  `Neural`           = c("S100B", "PLP1"))

## ---- inputs ------------------------------------------------------------------
md <- fread(file.path(CLUST, sprintf("%s_cell_metadata.csv", SAMPLE)))
setnames(md, c("V1", "spatial_snn_res.1", "umap_1", "umap_2"),
         c("cell", "cluster", "umap1", "umap2"))
md[, cluster := as.character(cluster)]
md[, y := -y]                                  # image y grows downward
lv <- as.character(sort(as.integer(unique(md$cluster))))
md[, cluster := factor(cluster, levels = lv)]
n_of <- md[, .N, by = cluster][order(cluster)]
log("cells: ", nrow(md), "  res", RES, " clusters: ", length(lv))

asg <- fread(st("major_assignment.csv"))
asg[, cluster := as.character(cluster)]
stopifnot(setequal(asg$cluster, lv))
# a fixed, readable order: epithelium, then lymphoid, myeloid, mesenchymal, vascular
MAJ_ORDER <- c("Tumour / luminal", "Basal / myoepithelial", "Plasma cell", "T / NK",
               "Mast cell", "Myeloid", "Fibroblast / CAF", "Endothelial",
               "Pericyte / SMC", "Adipocyte")
MAJ_LV <- MAJ_ORDER[MAJ_ORDER %in% asg$major]
stopifnot(setequal(MAJ_LV, unique(asg$major)))
md[asg, major := i.major, on = "cluster"]
md[, major := factor(major, levels = MAJ_LV)]
log("major types: ", paste(MAJ_LV, collapse = ", "))

b  <- readRDS(st("markers_input.rds"))
bc <- factor(b$idents, levels = lv)
top <- fread(st("top10_markers.csv"))[, cluster := as.character(cluster)]

## ---- palettes -----------------------------------------------------------------
# 31 hues that stay apart on white at a 0.1 pt dot; the near-duplicates are spread
# so that clusters sharing a hue are not neighbours on the tissue.
PAL31 <- c("#E6194B","#3CB44B","#4363D8","#F58231","#911EB4","#008080","#9A6324",
           "#F032E6","#000075","#808000","#00CED1","#800000","#AAFFC3","#FFD8B1",
           "#A9A9A9","#FF4500","#1E90FF","#228B22","#DA70D6","#8B4513","#20B2AA",
           "#DC143C","#4682B4","#9ACD32","#FF1493","#2F4F4F","#CD853F","#6A5ACD",
           "#00FA9A","#B22222","#5F9EA0")
PAL <- setNames(PAL31[seq_along(lv)], lv)
MAJ_PAL <- setNames(c("#C44E52","#E377C2","#4C72B0","#17BECF","#9467BD","#8C564B",
                      "#BCBD22","#55A868","#2F6F9F","#7F7F7F","#FF9896","#AEC7E8",
                      "#D62728")[seq_along(MAJ_LV)], MAJ_LV)

th_sp <- theme_void(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 8, face = "bold", hjust = .5, margin = margin(b = 2)),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"))
asp <- diff(range(md$y)) / diff(range(md$x))

png_dev <- function(f, w, h) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(f, width = w, height = h, units = "in", res = 400, background = "white")
  else png(f, width = w, height = h, units = "in", res = 400, bg = "white", type = "cairo")
}

## ---- 2. all res1 clusters on the tissue ---------------------------------------
lab <- setNames(sprintf("%s  (%s)", lv,
                        formatC(n_of$N, big.mark = ",", format = "d")), lv)
# Every point the same size. Only the draw order helps the small clusters: they go
# last so they are not painted over. Sizing by rarity, as an earlier version did,
# makes a rare group look more abundant than it is.
md[, rare := frank(n_of$N[match(cluster, n_of$cluster)], ties.method = "dense")]
setorder(md, -rare)                            # small clusters drawn last, on top
p_cl <- ggplot(md, aes(x, y, colour = cluster)) +
  geom_point(size = PT, stroke = 0) +
  scale_colour_manual(values = PAL, labels = lab, name = sprintf("Cluster (res %s)", RES)) +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2.2), ncol = 2)) +
  labs(title = sprintf("%s res %s -- all %d clusters", SAMPLE, RES, length(lv))) + th_sp
W <- 10; H <- (W - 2.6) * asp + 0.5
png_dev(st("clusters_spatial.png"), W, H); print(p_cl); invisible(dev.off())
log("wrote clusters_spatial")

bg <- md[, .(x, y)]
p_fac <- wrap_plots(lapply(lv, function(g) {
  ggplot() +
    geom_point(data = bg, aes(x, y), colour = "grey90", size = 0.03, stroke = 0) +
    geom_point(data = md[cluster == g], aes(x, y), colour = PAL[[g]], size = PT, stroke = 0) +
    coord_fixed() + labs(title = lab[[g]]) + th_sp +
    theme(plot.title = element_text(size = 6.5, face = "bold", hjust = .5, colour = PAL[[g]]))
}), ncol = 6)
nr <- ceiling(length(lv) / 6); W2 <- 11; H2 <- nr * ((W2 / 6) * asp + 0.2)
png_dev(st("clusters_facets.png"), W2, H2); print(p_fac); invisible(dev.off())
log("wrote clusters_facets")

## ---- dot plot helper -----------------------------------------------------------
dotplot <- function(genes, grp, grp_lv, gene_split = NULL, title, sub) {
  genes <- genes[genes %in% rownames(b$data)]
  X <- b$data[genes, , drop = FALSE]
  pct <- t(sapply(grp_lv, function(g) Matrix::rowMeans(X[, grp == g, drop = FALSE] > 0))) * 100
  avg <- t(sapply(grp_lv, function(g) Matrix::rowMeans(X[, grp == g, drop = FALSE])))
  Z   <- scale(avg); Z[!is.finite(Z)] <- 0     # z-score each gene across groups
  D <- as.data.table(reshape2::melt(Z, varnames = c("group", "gene"), value.name = "z"))
  D[, pct := as.vector(reshape2::melt(pct)$value)]
  D[, gene := factor(gene, levels = genes)]
  D[, group := factor(group, levels = rev(grp_lv))]
  if (!is.null(gene_split)) D[, block := factor(gene_split[as.character(gene)],
                                                levels = unique(gene_split))]
  p <- ggplot(D, aes(gene, group, size = pct, colour = z)) +
    geom_point(stroke = 0) +
    scale_size_continuous(range = c(0.1, 3), limits = c(0, 100),
                          breaks = c(10, 50, 90), name = "% expressing") +
    scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0,
                           name = "Mean expr.\n(z across groups)") +
    labs(x = NULL, y = NULL, title = title, subtitle = sub) +
    theme_bw(base_size = 7) +
    theme(panel.grid.major = element_line(linewidth = .15, colour = "grey92"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 5),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.5, colour = "grey30"),
          legend.key.size = unit(.3, "cm"), legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6),
          strip.text = element_text(size = 5.5, margin = margin(1, 1, 1, 1)),
          plot.background = element_rect(fill = "white", colour = NA))
  if (!is.null(gene_split)) p <- p + facet_grid(cols = vars(block), scales = "free_x",
                                                space = "free_x") +
    theme(panel.spacing.x = unit(1.2, "pt"))
  p
}

## ---- 3. top markers x res1 cluster ---------------------------------------------
tt <- top[rank <= TOP_N][order(as.integer(cluster), rank)]
tt <- tt[!duplicated(gene)]
p_dot <- dotplot(tt$gene, bc, lv, setNames(tt$cluster, tt$gene),
                 sprintf("%s res %s -- top %d markers per cluster", SAMPLE, RES, TOP_N),
                 sprintf("%d unique genes; column blocks are the cluster each gene marks; %d cells per cluster",
                         nrow(tt), min(table(bc))))
ggsave(st("top_marker_dotplot.pdf"), p_dot, width = 0.085 * nrow(tt) + 2.4,
       height = 0.17 * length(lv) + 1.6, limitsize = FALSE)
ggsave(st("top_marker_dotplot.png"), p_dot, width = 0.085 * nrow(tt) + 2.4,
       height = 0.17 * length(lv) + 1.6, dpi = 300, limitsize = FALSE)
log("wrote top_marker_dotplot")

## ---- 4+5: tissue, UMAP and classic-marker dot plot, for one grouping ----------
# Run twice: once over the ten predicted types, once over the four compartments they
# collapse into. Rare groups are drawn last and larger, or 302 pericytes vanish under
# 450,941 tumour cells.
cg  <- unlist(CLASSIC, use.names = FALSE)
blk <- setNames(rep(names(CLASSIC), lengths(CLASSIC)), cg)
cg_present <- cg[cg %in% rownames(b$data)]
miss <- setdiff(cg, cg_present)
if (length(miss)) log("classic markers absent from the matrix: ", paste(miss, collapse = ", "))

figs_for <- function(col, grp_lv, pal, stem, what) {
  # Every point the same size: a dot's area must mean one cell and nothing else.
  # Only the draw order is tuned -- ggplot paints in row order, so the rare groups go
  # last and are not painted over by the abundant ones.
  sz <- md[, .N, by = c(col)][order(N)]
  md[, draw := match(as.character(get(col)), as.character(sz[[col]]))]
  setorder(md, -draw)
  cnt <- md[, .N, by = c(col)]
  lb  <- setNames(sprintf("%s  (%s)", cnt[[col]], formatC(cnt$N, big.mark = ",", format = "d")),
                  as.character(cnt[[col]]))
  p_s <- ggplot(md, aes(x, y, colour = .data[[col]])) +
    geom_point(size = PT, stroke = 0) +
    scale_colour_manual(values = pal, labels = lb[grp_lv], breaks = grp_lv, name = what) +
    coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 2.4), ncol = 1)) +
    labs(title = paste(what, "-- tissue")) + th_sp
  p_u <- ggplot(md, aes(umap1, umap2, colour = .data[[col]])) +
    geom_point(size = PT, stroke = 0, show.legend = FALSE) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(title = paste(what, "-- UMAP"), x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = .5),
          plot.background = element_rect(fill = "white", colour = NA))
  W3 <- 11; H3 <- max((W3 * 0.55 - 1.8) * asp, 4.2) + 0.4
  png_dev(st(paste0(stem, "_spatial_umap.png")), W3, H3)
  print(p_s + p_u + plot_layout(widths = c(1.25, 1))); invisible(dev.off())
  log("wrote ", stem, "_spatial_umap")

  bg2 <- factor(asg[[col]][match(as.character(bc), asg$cluster)], levels = grp_lv)
  p_c <- dotplot(cg_present, bg2, grp_lv, blk[cg_present],
                 sprintf("%s -- canonical lineage markers by %s", SAMPLE, tolower(what)),
                 "genes grouped by the lineage they mark; these are NOT the genes the calls were made from")
  w <- 0.1 * length(cg_present) + 2.6; h <- 0.22 * length(grp_lv) + 2.0
  ggsave(st(paste0(stem, "_classic_dotplot.pdf")), p_c, width = w, height = h, limitsize = FALSE)
  ggsave(st(paste0(stem, "_classic_dotplot.png")), p_c, width = w, height = h, dpi = 300, limitsize = FALSE)
  log("wrote ", stem, "_classic_dotplot")
}

figs_for("major", MAJ_LV, MAJ_PAL, "major", "Predicted major cell type")

## ---- the four compartments the ten types collapse into -------------------------
COMPARTMENT <- c("Tumour / luminal"      = "Epithelial",
                 "Basal / myoepithelial" = "Epithelial",
                 "Plasma cell"           = "Immune",
                 "T / NK"                = "Immune",
                 "Mast cell"             = "Immune",
                 "Myeloid"               = "Immune",
                 "Fibroblast / CAF"      = "Stroma / SMC",
                 "Endothelial"           = "Stroma / SMC",
                 "Pericyte / SMC"        = "Stroma / SMC",
                 "Adipocyte"             = "Adipocyte")
COMP_LV  <- c("Epithelial", "Immune", "Stroma / SMC", "Adipocyte")
COMP_PAL <- c(Epithelial = "#C44E52", Immune = "#4C72B0",
              `Stroma / SMC` = "#55A868", Adipocyte = "#DD8452")
asg[, compartment := COMPARTMENT[major]]
md[, compartment := factor(COMPARTMENT[as.character(major)], levels = COMP_LV)]
stopifnot(!anyNA(md$compartment))
print(md[, .N, by = compartment][order(match(compartment, COMP_LV))])
figs_for("compartment", COMP_LV, COMP_PAL, "compartment", "Compartment")

## ---- epithelium against everything else ---------------------------------------
# The same three figures once more at the coarsest grouping, which is the one the
# CNV work actually needs: inferCNV runs on the epithelium and references the rest.
EPI_LV  <- c("Epithelium", "Rest")
# Rest is 36% of the section, not background, so it gets a real colour rather
# than the grey an "other" category usually gets.
EPI_PAL <- c(Epithelium = "#C44E52", Rest = "#4C72B0")
asg[, epi := fifelse(compartment == "Epithelial", "Epithelium", "Rest")]
md[, epi := factor(fifelse(compartment == "Epithelial", "Epithelium", "Rest"), levels = EPI_LV)]
print(md[, .N, by = epi][order(match(epi, EPI_LV))])
figs_for("epi", EPI_LV, EPI_PAL, "epi_vs_rest", "Epithelium vs rest")

fwrite(md[, .(cell, x, y, umap1, umap2, cluster, major, compartment, epi)],
       st("cells_with_major.csv"))
fwrite(asg, st("major_assignment_with_compartment.csv"))
log("done -> ", OUT)
