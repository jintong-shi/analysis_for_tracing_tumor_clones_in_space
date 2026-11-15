#!/usr/bin/env Rscript
# =====================================================================
# Fig 5 (high-resolution ST) — MAIN composed figure, rebuilt from scratch
# using ONLY high_res_st/input_data. Every panel is replotted from data;
# every plotted value (schematic geometry, pseudobulk means, correlations,
# subclone pseudobulk, all-epi binning) is recomputed here. Sources ONLY
# utils/high_res_setup.R.
#
# Panels (letters = reading order). Rows 3-4 pair one CNV heatmap with the matching
# spatial map, once per platform: Xenium (d/e), StrataMap (f/g). Row 5 is one 3x3
# tool-concordance square per platform.
#   a schematic (2/8/16um bins, segmented cells, meta-cells)      [schematic/]
#   b pseudobulk Pearson r vs GT bar                              [visiumHD/infercnv]
#   c subclone correspondence (HD r40 k2 vs Visium k3)            [visiumHD/subclone + infercnv]
#   d atera benign + subclone CNV heatmap                         [atera/subclone_cnv_input]
#   e atera normal + subclone spatial                             [atera/cells + labels]
#   f StrataMap G3 reference + CNV-cluster heatmap                [novaseq/cnv_input]
#   g StrataMap G3 normal + CNV-cluster spatial                   [novaseq/coords_labels]
#   h/i/j inferCNV-CopyKAT-CONICSmat concordance, per platform    [tool_concordance/]
#
# The 4-strategy pseudobulk inferCNV heatmap that used to be panel b now lives in
# the VisiumHD supplement (supp_visiumHD.R panel j); the correlation bar that
# summarises it stayed here as panel b.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     paper/figures_final/high_res_st/scripts/fig5_main.R
# =====================================================================
suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(cowplot)
  library(ComplexHeatmap); library(circlize); library(grid)
}))

## ---- locate self / paths --------------------------------------------
.this <- (function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) normalizePath(sub("^--file=", "", m[1])) else normalizePath("fig5_main.R")
})()
SDIR  <- dirname(.this)
TDIR  <- normalizePath(file.path(SDIR, ".."))
UTILS <- normalizePath(file.path(TDIR, "..", "utils"))
IN    <- file.path(TDIR, "input_data")
OUTF  <- file.path(TDIR, "main_figures"); dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)
OUTT  <- file.path(TDIR, "data_tables");  dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS, "high_res_setup.R"))
wtab <- function(df, nm) fwrite(as.data.frame(df), file.path(OUTT, nm))

FONT <- pick_font()
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))
message("font: ", FONT)

## ---- shared style (FIGURE_STANDARDS 0.1) ----------------------------
S_TIT <- 6.5; S_AXT <- 6; S_TICK <- 5.5; S_LEGT <- 6; S_LET <- 8
CNV_COL <- cnv_colorramp()
theme_f5 <- function(void = FALSE) {
  base <- if (void) theme_void(base_size = S_TICK, base_family = FONT)
          else theme_classic(base_size = S_TICK, base_family = FONT)
  base + theme(text = element_text(family = FONT),
    plot.title = element_text(family = FONT, face = "bold", size = S_TIT, hjust = 0.5, margin = margin(b = 1)),
    axis.title = element_text(size = S_AXT), axis.title.x = element_text(margin = margin(t = 1)),
    axis.title.y = element_text(margin = margin(r = 1)),
    axis.text = element_text(size = S_TICK, colour = "black"),
    legend.title = element_text(size = S_LEGT), legend.text = element_text(size = S_TICK),
    legend.key.size = unit(0.24, "cm"), plot.margin = margin(2, 2, 2, 2))
}
heat_grob <- function(expr) grid.grabExpr(expr, warn = 0)

## save one panel standalone (handles ggplot OR ComplexHeatmap grob) --------
save_panel <- function(p, stem, w, h) save_panel_both(p, stem, OUTF, w, h, FONT)

# ============================================================ a: schematic
build_a <- function() {
  b  <- fread(file.path(IN, "schematic", "schematic_window_bins.csv"))
  bd <- fread(file.path(IN, "schematic", "schematic_window_bounds.csv"))
  xr <- c(bd$xmin, bd$xmax); yr <- c(bd$ymin, bd$ymax)
  b[, is_epi := as.logical(is_epi)]
  b[, celltype := ifelse(is_epi, "Epithelial", "Other")]
  EPI <- c("Epithelial" = "#4DAF4A", "Other" = "grey80"); grey_col <- "grey80"

  ## 8um bins (majority-rule epithelial, as in the figure) -------------
  bin8 <- b[, .(n = .N, ne = sum(is_epi & oid != 0), nl = sum(oid != 0)),
            by = .(x8 = x %/% 4, y8 = y %/% 4)][n >= 8]
  bin8[, `:=`(cx = 4 * x8 + 1.5, cy = 4 * y8 + 1.5,
              celltype = ifelse(nl > 0 & ne / nl >= 0.5, "Epithelial", "Other"))]
  ## 16um bins (epithelial fraction) ----------------------------------
  bin16 <- b[, .(n = .N, ne = sum(is_epi & oid != 0), nl = sum(oid != 0)),
             by = .(x16 = x %/% 8, y16 = y %/% 8)][n >= 32]
  bin16[, `:=`(cx = 8 * x16 + 3.5, cy = 8 * y16 + 3.5, frac = fifelse(nl > 0, ne / nl, 0))]
  ## cell outline segments --------------------------------------------
  cc <- b[oid != 0, .(x, y, cid = oid)]; setkey(cc, x, y)
  nb <- function(dx, dy) cc[.(cc$x + dx, cc$y + dy), cid, on = .(x, y)]
  segs <- rbind(
    cc[cid != nb(1, 0)  | is.na(nb(1, 0)),  .(x = x + .5, xend = x + .5, y = y - .5, yend = y + .5)],
    cc[cid != nb(-1, 0) | is.na(nb(-1, 0)), .(x = x - .5, xend = x - .5, y = y - .5, yend = y + .5)],
    cc[cid != nb(0, 1)  | is.na(nb(0, 1)),  .(x = x - .5, xend = x + .5, y = y + .5, yend = y + .5)],
    cc[cid != nb(0, -1) | is.na(nb(0, -1)), .(x = x - .5, xend = x + .5, y = y - .5, yend = y - .5)])
  ## meta-cell membership + cleanup -----------------------------------
  wlab <- b[oid != 0]
  wlab[, meta_show := NA_character_]
  wlab[is_epi & !is.na(pseudo_id), meta_show := as.character(pseudo_id)]
  largest_component <- function(cx, cy, thr = 12) {
    n <- length(cx); if (n == 1) return(1L)
    adj <- as.matrix(dist(cbind(cx, cy))) <= thr; comp <- integer(n); k <- 0L
    for (i in seq_len(n)) if (comp[i] == 0L) { k <- k + 1L; st <- i
      while (length(st)) { j <- st[1]; st <- st[-1]
        if (comp[j] == 0L) { comp[j] <- k; st <- c(st, which(adj[j, ] & comp == 0L)) } } }
    comp }
  cellc <- wlab[!is.na(meta_show), .(cx = mean(x), cy = mean(y)), by = .(meta_show, oid)]
  cellc[, comp := largest_component(cx, cy), by = meta_show]
  keep <- cellc[, .SD[comp == as.integer(names(which.max(table(comp))))], by = meta_show]$oid
  wlab[!is.na(meta_show) & !(oid %in% keep), meta_show := NA_character_]
  small <- wlab[!is.na(meta_show), uniqueN(oid), by = meta_show][V1 < 4]$meta_show
  wlab[meta_show %in% small, meta_show := NA_character_]
  hulls <- wlab[!is.na(meta_show), { h <- chull(x, y); .(hx = x[h], hy = y[h]) }, by = meta_show]
  hulls[, `:=`(mx = mean(hx), my = mean(hy)), by = meta_show]
  dd <- sqrt((hulls$hx - hulls$mx)^2 + (hulls$hy - hulls$my)^2) + 1e-9
  hulls[, `:=`(hx = hx + 1.2 * (hx - mx) / dd, hy = hy + 1.2 * (hy - my) / dd)]
  metas_in <- sort(unique(na.omit(wlab$meta_show)))
  qual <- setdiff(unique(c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(9, "Set1"),
                           RColorBrewer::brewer.pal(8, "Set2"))), c("#4DAF4A", "#999999", "#666666", "#B3B3B3"))
  set.seed(2); pal <- setNames(sample(qual, length(metas_in)), metas_in)

  ## data tables ------------------------------------------------------
  wtab(b[, .(x, y, oid, is_epi, pseudo_id, celltype)], "fig5_panel_a_bins2um.csv")
  wtab(bin8[, .(cx, cy, ne, nl, celltype)], "fig5_panel_a_bins8um.csv")
  wtab(bin16[, .(cx, cy, ne, nl, frac)], "fig5_panel_a_bins16um.csv")
  wtab(wlab[!is.na(meta_show), .(x, y, oid, meta_show)], "fig5_panel_a_metacell_bins.csv")

  th <- theme_void(base_size = 5, base_family = FONT) +
    theme(plot.title = element_text(hjust = .5, face = "bold", size = 5.5, family = FONT),
          legend.position = "none", plot.margin = margin(1, 1, 1, 1))
  lims <- list(coord_equal(xlim = xr, ylim = yr, expand = FALSE))
  pA <- ggplot(b, aes(x, y)) + geom_tile(width = .85, height = .85, fill = "grey80") + lims + th + ggtitle("2 µm bins (raw)")
  pB <- ggplot(bin8) + geom_tile(aes(cx, cy, fill = celltype), width = 4, height = 4, colour = "grey30", linewidth = .2) +
    scale_fill_manual(name = NULL, values = EPI) + lims + th + ggtitle("8 µm bins")
  pB16 <- ggplot(bin16) + geom_tile(aes(cx, cy, fill = frac), width = 8, height = 8, colour = "grey30", linewidth = .2) +
    scale_fill_gradient(name = "Epi frac", low = "#EDEDED", high = "#1B5E20", limits = c(0, 1), breaks = c(0, .5, 1)) +
    lims + th + ggtitle("16 µm bins")
  pC <- ggplot() +
    geom_tile(data = b[oid == 0], aes(x, y), width = 1, height = 1, fill = grey_col) +
    geom_tile(data = b[oid != 0], aes(x, y, fill = celltype), width = 1, height = 1) +
    geom_segment(data = segs, aes(x, xend = xend, y = y, yend = yend), colour = "white", linewidth = .25) +
    scale_fill_manual(name = NULL, values = EPI) + lims + th + ggtitle("Segmented cells")
  pD <- ggplot() +
    geom_tile(data = b[oid == 0], aes(x, y), width = 1, height = 1, fill = grey_col) +
    geom_tile(data = wlab[is.na(meta_show)], aes(x, y), width = 1, height = 1, fill = grey_col) +
    geom_tile(data = wlab[!is.na(meta_show)], aes(x, y, fill = meta_show), width = 1, height = 1) +
    geom_segment(data = segs, aes(x, xend = xend, y = y, yend = yend), colour = "white", linewidth = .22) +
    geom_polygon(data = hulls, aes(hx, hy, group = meta_show, colour = meta_show), fill = NA, linewidth = .5) +
    scale_fill_manual(values = pal, na.value = NA) + scale_colour_manual(values = pal, guide = "none") +
    lims + th + ggtitle("Meta-cells (r40)")
  squares <- plot_grid(pA, pB, pB16, pC, pD, nrow = 1, ncol = 5)
  leg_epi <- cowplot::get_legend(
    ggplot(data.frame(x = 1:2, ct = factor(c("Epithelial", "Other"), levels = c("Epithelial", "Other")))) +
      geom_tile(aes(x, 1, fill = ct)) + scale_fill_manual(name = NULL, values = EPI) +
      theme(legend.position = "bottom", legend.text = element_text(size = 5.5, family = FONT),
            legend.key.size = unit(0.20, "cm"), legend.margin = margin(0, 0, 0, 0)))
  leg_frac <- cowplot::get_legend(
    ggplot(data.frame(x = 1:2, frac = c(0, 1))) + geom_tile(aes(x, 1, fill = frac)) +
      scale_fill_gradient(name = "Epi frac", low = "#EDEDED", high = "#1B5E20", limits = c(0, 1), breaks = c(0, .5, 1)) +
      theme(legend.position = "bottom", legend.title = element_text(size = 5.5, family = FONT),
            legend.text = element_text(size = 5, family = FONT),
            legend.key.height = unit(0.16, "cm"), legend.key.width = unit(0.45, "cm"), legend.margin = margin(0, 0, 0, 0)))
  plot_grid(squares, plot_grid(leg_epi, leg_frac, nrow = 1), ncol = 1, rel_heights = c(1, 0.12))
}

# ============================================================ b: correlation bar
build_b <- function() {
  cc <- copy(compute_pseudobulk(IN)$corr)[order(pearson)]
  cc[, strategy := factor(strategy, levels = strategy)]
  wtab(cc[, .(strategy = LABEL_MAP[as.character(strategy)], n_shared, pearson)], "fig5_panel_b_pseudobulk_corr.csv")
  labs <- setNames(sprintf("%s\n(n=%d)", LABEL_MAP[as.character(cc$strategy)], cc$n_shared), cc$strategy)
  ggplot(cc, aes(strategy, pearson, fill = strategy)) +
    geom_col(width = .8, colour = "white") +
    geom_text(aes(label = sprintf("%.2f", pearson)), hjust = -0.2, size = 5.5 / .pt, fontface = "bold") +
    scale_fill_manual(values = HD_COLS, guide = "none") + scale_x_discrete(labels = labs) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, .5), expand = expansion(c(0, 0))) +
    coord_flip() + labs(x = NULL, y = "Pearson r vs GT") + theme_f5() +
    theme(axis.text.y = element_text(size = 5, lineheight = 0.85))
}

# ============================================================ c: subclone correspondence
build_c <- function() {
  SB <- file.path(IN, "visiumHD", "subclone")
  r <- readRDS(file.path(SB, "prep_r40_k2.rds"))
  v <- readRDS(file.path(SB, "prep_visium_k3.rds"))
  r_xy <- fread(file.path(SB, "coords_r40.csv")); setnames(r_xy, c("id", "x", "y"))
  v_xy <- fread(file.path(SB, "coords_visium.csv")); setnames(v_xy, c("id", "x", "y"))
  # benign (inferCNV "Normal" reference) spots for the Visium map background
  gt <- read_infercnv(file.path(IN, "visiumHD", "infercnv", "visium_pt10_AP5_GT_run.final.infercnv_obj"))
  v_benign <- colnames(gt@expr.data)[sort(unique(as.integer(unlist(gt@reference_grouped_cell_indices))))]
  # recompute subclone pseudobulk from the raw spot x gene matrices
  pb_by_sub <- function(mat, sub) {
    sub <- sub[rownames(mat)]; lv <- sort(unique(sub))
    t(sapply(lv, function(s) colMeans(mat[sub == s, , drop = FALSE]))) }
  rpb <- pb_by_sub(r$mat, r$sub); vpb <- pb_by_sub(v$mat, v$sub)
  r_sub <- r$sub; v_sub <- v$sub
  r40_map <- c(subclone_1 = "Clone 1", subclone_2 = "Clone 2")
  vis_map <- c(subclone_1 = "Clone 1", subclone_3 = "Clone 2", subclone_2 = "Clone 3")
  MAP_ASP <- 2.05
  smap <- function(xy, sub, cmap, title, pt, bg_ids = NULL, box_asp = MAP_ASP) {
    co <- as.data.frame(xy); co$Clone <- cmap[sub[co$id]]
    fg <- co[!is.na(co$Clone), ]; fg$Clone <- factor(fg$Clone, levels = names(CLONE))
    bg <- if (is.null(bg_ids)) co[is.na(co$Clone), ] else co[is.na(co$Clone) & co$id %in% bg_ids, ]
    shown <- rbind(bg[, c("x", "y")], fg[, c("x", "y")])
    L <- spatial_lims(shown$x, -shown$y, box_asp)
    ggplot() +
      geom_point(data = bg, aes(x, -y), color = "grey85", size = pt, stroke = 0) +
      geom_point(data = fg, aes(x, -y, color = Clone), size = pt, stroke = 0) +
      scale_color_manual(values = CLONE, drop = FALSE, guide = "none") +
      coord_fixed(xlim = L$x, ylim = L$y, expand = FALSE) + labs(title = title) +
      theme_void(base_size = 5, base_family = FONT) +
      theme(plot.title = element_text(face = "bold", hjust = .5, size = 5.5, family = FONT), plot.margin = margin(1, 1, 1, 1)) }
  m_r40 <- smap(r_xy, r_sub, r40_map, "HD epi-meta-cells (r40, k=2)", .5)
  m_vis <- smap(v_xy, v_sub, vis_map, "regular-Visium (pt10, k=3)", .7, bg_ids = v_benign)
  # data tables: spatial + pseudobulk
  dt_map <- function(xy, sub, cmap, ds) { co <- as.data.frame(xy); co$clone <- cmap[sub[co$id]]
    data.table(dataset = ds, id = co$id, x = co$x, y = co$y, clone = co$clone) }
  wtab(rbind(dt_map(r_xy, r_sub, r40_map, "HD epi-meta-cells"),
             dt_map(v_xy, v_sub, vis_map, "regular-Visium")), "fig5_panel_c_subclone_maps.csv")
  ## pseudobulk stack on common genes
  common <- colnames(r$mat)[colnames(r$mat) %in% colnames(vpb)]
  chr_c <- droplevels(factor(as.character(r$chr_genes), levels = levels(r$chr_genes))[match(common, colnames(r$mat))])
  n_sl <- length(levels(chr_c)); levels(chr_c) <- short_chr(levels(chr_c))
  rownames(rpb) <- r40_map[rownames(rpb)]; rpb <- rpb[c("Clone 1", "Clone 2"), common, drop = FALSE]
  rownames(vpb) <- vis_map[rownames(vpb)]; vpb <- vpb[c("Clone 1", "Clone 2", "Clone 3"), common, drop = FALSE]
  lr <- reshape2::melt(rpb, varnames = c("clone", "gene"), value.name = "cnv"); lr$dataset <- "HD epi-meta-cells"
  lv <- reshape2::melt(vpb, varnames = c("clone", "gene"), value.name = "cnv"); lv$dataset <- "regular-Visium"
  pball <- as.data.table(rbind(lr, lv)); pball[, chr := as.character(chr_c)[match(gene, common)]]
  wtab(pball[, .(dataset, clone, gene, chr, cnv)], "fig5_panel_c_subclone_pseudobulk.csv")
  pbht <- function(mat, nm, show_chr, show_leg) {
    a <- list(matrix = mat, name = nm, col = CNV_COL, cluster_columns = FALSE, cluster_rows = FALSE,
      row_order = rownames(mat), show_row_names = TRUE, row_names_side = "left",
      row_names_gp = gpar(fontsize = S_TICK, fontfamily = FONT), show_column_names = FALSE,
      column_split = chr_c, column_gap = unit(0, "mm"), column_title_gp = gpar(fontsize = 5, fontfamily = FONT),
      height = unit(nrow(mat), "null"), border = TRUE, show_heatmap_legend = show_leg,
      left_annotation = rowAnnotation(Clone = factor(rownames(mat), levels = names(CLONE)), col = list(Clone = CLONE),
        show_annotation_name = FALSE, show_legend = show_leg, simple_anno_size = unit(2, "mm"),
        annotation_legend_param = list(Clone = list(title = "Subclone", title_gp = gpar(fontsize = S_LEGT, fontfamily = FONT),
          labels_gp = gpar(fontsize = S_TICK, fontfamily = FONT), grid_height = unit(2.6, "mm")))),
      heatmap_legend_param = list(title = "inferCNV", at = CNV_BREAKS, title_gp = gpar(fontsize = S_LEGT, fontfamily = FONT),
        labels_gp = gpar(fontsize = S_TICK, fontfamily = FONT), grid_width = unit(2.2, "mm")))
    if (!show_chr) a$column_title <- rep("", n_sl)
    do.call(Heatmap, a) }
  stack <- pbht(rpb, "inferCNV", TRUE, TRUE) %v% pbht(vpb, "inferCNV2", FALSE, FALSE)
  g_ht <- heat_grob(draw(stack, ht_gap = unit(2, "mm"), merge_legend = TRUE, heatmap_legend_side = "right",
                         padding = unit(c(1, 1, 1, 1), "mm")))
  ggdraw() + draw_plot(m_r40, 0.00, 0.50, 0.26, 0.50) + draw_plot(m_vis, 0.00, 0.00, 0.26, 0.50) +
    draw_plot(g_ht, 0.27, 0.00, 0.73, 1.00)
}

# ================================================ shared: clone colours + panels
## One palette per platform; the two are not matched, so nothing suggests that a Xenium
## clone corresponds to a StrataMap cluster. Normal/benign is green in both. Within Xenium
## the four hues are spread so the adjacent pairs A/C and B/D stay separable at point size.
NORMAL_COL      <- "#55A868"                                  # benign / reference cells
PT              <- 0.12                                       # spatial point size
CLONE_COLS      <- c("#C1521D", "#7D3C98", "#3C8DBC", "#E8B21A")   # Xenium A-D
CLONE_COLS_NOVA <- c("#2C5F9E", "#C2185B")                         # StrataMap C1-C2
pal_for <- function(clones, cols = CLONE_COLS) setNames(cols[seq_along(clones)], clones)

## reference/normal cells as the top block, observation rows below split by clone.
cnv_panel <- function(mat_ref, mat_obs, obs_lab, chr_genes, ref_lab, clones,
                      cols = CLONE_COLS) {
  cg     <- factor(short_chr(as.character(chr_genes)),
                   levels = unique(short_chr(as.character(chr_genes))))
  grp_lv <- c(ref_lab, clones)
  grp_c  <- c(setNames(NORMAL_COL, ref_lab), pal_for(clones, cols))
  common <- function(...) list(col = CNV_COL, cluster_columns = FALSE, show_row_names = FALSE,
    show_column_names = FALSE, column_split = cg, column_gap = unit(0, "mm"),
    column_title_gp = gpar(fontsize = 5, fontfamily = FONT), border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.3),
    use_raster = TRUE, raster_quality = 2, row_title_rot = 0,
    row_title_gp = gpar(fontsize = S_AXT, fontface = "bold", fontfamily = FONT), ...)
  ht_r <- do.call(Heatmap, c(list(mat_ref, name = "inferCNV", row_title = NULL,
    cluster_rows = FALSE, height = unit(1, "null"),
    left_annotation = rowAnnotation(
      Group = factor(rep(ref_lab, nrow(mat_ref)), levels = grp_lv), col = list(Group = grp_c),
      show_annotation_name = FALSE, simple_anno_size = unit(2, "mm"),
      annotation_legend_param = list(Group = list(title = "Group", ncol = 1, direction = "vertical",
        title_gp = gpar(fontsize = S_LEGT, fontfamily = FONT),
        labels_gp = gpar(fontsize = S_TICK, fontfamily = FONT),
        grid_height = unit(2.6, "mm"), grid_width = unit(2.6, "mm")))),
    heatmap_legend_param = list(title = "inferCNV", at = CNV_BREAKS, direction = "vertical",
      title_gp = gpar(fontsize = S_LEGT, fontfamily = FONT),
      labels_gp = gpar(fontsize = S_TICK, fontfamily = FONT),
      grid_width = unit(2.2, "mm"), legend_height = unit(12, "mm"))), common()))
  ht_o <- do.call(Heatmap, c(list(mat_obs, name = "obs", show_heatmap_legend = FALSE,
    row_title = NULL, row_split = factor(obs_lab, levels = clones),
    cluster_rows = FALSE, cluster_row_slices = FALSE,
    height = unit(5, "null"), row_gap = unit(0, "mm"),
    left_annotation = rowAnnotation(Group = factor(obs_lab, levels = grp_lv),
      col = list(Group = grp_c), show_annotation_name = FALSE, show_legend = FALSE,
      simple_anno_size = unit(2, "mm"))), common()))
  heat_grob(draw(ht_r %v% ht_o, ht_gap = unit(2, "mm"), merge_legend = TRUE,
                 heatmap_legend_side = "right", annotation_legend_side = "right",
                 padding = unit(c(1, 1, 1, 1), "mm")))
}

## spatial map: normal cells plus one colour per clone. Normal is far rarer than the
## clones on both platforms, so it is a separate layer drawn last (on top); ref_size
## enlarges it where the count is too small to read at the clone point size (StrataMap
## has 476 normal cells against 90,000 tumour ones). The legend keys stay uniform.
spatial_panel <- function(df, clones, box_asp, ref_lab, legend_title,
                          cols = CLONE_COLS, ref_size = PT) {
  df[, group := factor(group, levels = c(ref_lab, clones))]
  L <- spatial_lims(df$x, df$y, box_asp)
  ggplot(mapping = aes(x, y, colour = group)) +
    geom_point(data = df[group != ref_lab], size = PT, stroke = 0) +
    geom_point(data = df[group == ref_lab], size = ref_size, stroke = 0) +
    scale_colour_manual(values = c(setNames(NORMAL_COL, ref_lab), pal_for(clones, cols)),
                        name = legend_title, drop = FALSE) +
    scale_x_continuous(limits = L$x, expand = c(0, 0)) +
    scale_y_continuous(limits = L$y, expand = c(0, 0)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8))) +
    theme_f5(void = TRUE) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = S_LEGT), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(1, 1, 1, 1),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
}

# ============================================================ d: atera CNV heatmap
build_d <- function() {
  d <- readRDS(file.path(IN, "atera", "atera_subclone_cnv_input.rds"))
  clones <- sort(unique(d$sub_epi))
  BIN <- 34L
  ## average every BIN consecutive rows of a clone, so the rasteriser sees a few
  ## thousand rows instead of ~76,000. Rows were shuffled in the prep, so a bin is
  ## an unbiased mean of that clone, not its most similar cells.
  bin_within <- function(mat, labels, lv, k) {
    parts <- lapply(lv, function(l) { mm <- mat[labels == l, , drop = FALSE]
      g <- ceiling(seq_len(nrow(mm)) / k); rowsum(mm, g) / as.numeric(table(g)) })
    labs <- unlist(lapply(lv, function(l) rep(l, ceiling(sum(labels == l) / k))))
    list(mat = do.call(rbind, parts), lab = labs) }
  epi <- bin_within(d$mat_epi, d$sub_epi, clones, BIN)
  wtab(data.table(gene = colnames(d$mat_epi), chr = as.character(d$chr_genes)), "fig5_panel_d_gene_annot.csv")
  wtab(data.table(row_id = seq_len(nrow(epi$mat)), group = epi$lab), "fig5_panel_d_epi_row_annot.csv")
  fwrite(as.data.table(round(d$mat_benign, 4)), file.path(OUTT, "fig5_panel_d_benign_matrix.csv.gz"))
  fwrite(as.data.table(round(epi$mat, 4)), file.path(OUTT, "fig5_panel_d_epi_binned_matrix.csv.gz"))
  message(sprintf("panel d: benign=%d epi=%d -> %d binned rows, clones %s",
                  d$n_benign, d$n_epi, nrow(epi$mat), paste(clones, collapse = ",")))
  cnv_panel(d$mat_benign, epi$mat, epi$lab, d$chr_genes, "Benign", clones)
}

# ============================================================ e: atera spatial
build_e <- function(box_asp) {
  co  <- fread(file.path(IN, "atera", "cells_coords.csv"))
  hl  <- fread(file.path(IN, "atera", "epi_cancer_normal_labels.csv"))
  cl  <- fread(file.path(IN, "atera", "atera_subclone_cell_labels.csv"))
  df  <- merge(co, cl, by = "cell_id", all.x = TRUE)
  df  <- merge(df, hl, by = "cell_id", all.x = TRUE)
  df[, group := ifelse(!is.na(clone), clone, ifelse(label == "Normal", "Normal", NA_character_))]
  df <- df[!is.na(group)]
  clones <- sort(unique(cl$clone))
  setnames(df, c("x_centroid", "y_centroid"), c("x", "y"))
  df[, y := -y]
  wtab(df[, .(cell_id, x, y, group)], "fig5_panel_e_atera_spatial.csv")
  message("panel e: ", paste(sprintf("%s=%d", names(table(df$group)), table(df$group)), collapse = " "))
  spatial_panel(df, clones, box_asp, "Normal", "Xenium clone")
}

# StrataMap C3 folded into C2 (Jintong 2026-10-02): the two were not distinct enough to
# keep apart, matching the subclone merge in high_res_ST/code/cnv_bin_utils.R that
# feeds panel j. Relabelled here rather than in the prep so the frozen extracts stay
# as captured; the prep's 12-in-cluster bins never straddled a cluster, so merging the
# labels needs no re-binning, and ComplexHeatmap's row_split re-groups the rows.
nova_merge <- function(x) ifelse(x == "C3", "C2", x)

# ========================================================== f: StrataMap CNV heatmap
build_f <- function() {
  d <- readRDS(file.path(IN, "novaseq", "novaseq_g3_merged_m3_cnv_input.rds"))
  d$bin_lab <- nova_merge(d$bin_lab)           # obs rows are already binned by the prep
  clones <- sort(unique(d$bin_lab))
  # The merge leaves the old C3 rows stacked at the bottom of C2, where they read as a
  # band. Shuffle within C2 so the block shows the clone's spread rather than the order
  # the two halves arrived in. Rows are unclustered here (cluster_rows = FALSE), so row
  # order carries no information and only C2 moves; seeded so the figure is stable.
  set.seed(7); i2 <- which(d$bin_lab == "C2")
  d$mat_bin[i2, ] <- d$mat_bin[sample(i2), , drop = FALSE]
  wtab(data.table(gene = colnames(d$mat_bin), chr = as.character(d$chr_genes)), "fig5_panel_f_gene_annot.csv")
  wtab(data.table(row_id = seq_len(nrow(d$mat_bin)), group = d$bin_lab), "fig5_panel_f_obs_row_annot.csv")
  fwrite(as.data.table(round(d$mat_ref, 4)), file.path(OUTT, "fig5_panel_f_ref_matrix.csv.gz"))
  fwrite(as.data.table(round(d$mat_bin, 4)), file.path(OUTT, "fig5_panel_f_obs_binned_matrix.csv.gz"))
  message(sprintf("panel f: ref=%d obs=%d -> %d binned rows (bin=%d), clusters %s",
                  d$n_ref, d$n_obs, nrow(d$mat_bin), d$bin, paste(clones, collapse = ",")))
  cnv_panel(d$mat_ref, d$mat_bin, d$bin_lab, d$chr_genes, "Normal", clones, CLONE_COLS_NOVA)
}

# ============================================================ g: StrataMap spatial
build_g <- function(box_asp) {
  df <- fread(file.path(IN, "novaseq", "novaseq_g3_coords_labels.csv"))
  df[, group := nova_merge(group)]
  clones <- sort(setdiff(unique(df$group), "Normal"))
  df[, y := -y]
  wtab(df[, .(barcode, x, y, group)], "fig5_panel_g_strataMap_spatial.csv")
  message("panel g: ", paste(sprintf("%s=%d", names(table(df$group)), table(df$group)), collapse = " "))
  spatial_panel(df, clones, box_asp, "Normal", "StrataMap CNV cluster",
                CLONE_COLS_NOVA, ref_size = 0.6)
}

# ======================== h/i/j: CNV-tool concordance, one matrix per platform
# Each subclone is pseudobulked per tool (high_res_ST/code/subclone_tool_concordance.R,
# 10Mb bins, all three tools on the common log2 scale), the three tools are correlated
# within each subclone, and the 3x3 is averaged over that platform's subclones. So the
# subclones act as replicates here rather than as separate panels -- the question is
# whether the tools agree on a platform, not on any one subclone.
TOOL_PLATFORMS <- c(h = "visiumhd", i = "atera_ds20k", j = "novaseq_ds20k")
# Canonical dataset labels (paper/CLAUDE.md banlist: never "atera" in a figure label).
TOOL_TITLES    <- c(h = "Prostate (VisiumHD)", i = "Breast (Xenium)",
                    j = "Breast (StrataMap, IDC grade 3)")
build_tool_r <- function(key) {
  tc <- tool_concordance(IN, TOOL_PLATFORMS[[key]])
  d  <- as.data.table(reshape2::melt(tc$r, varnames = c("tool_a", "tool_b"), value.name = "r"))
  wtab(d, sprintf("fig5_panel_%s_tool_concordance_%s.csv", key, TOOL_PLATFORMS[[key]]))
  message(sprintf("panel %s: %s, %d subclones, %d bins", key, TOOL_PLATFORMS[[key]],
                  tc$n_sub, tc$n_bin))
  d[, `:=`(tool_a = factor(tool_a, levels = CNV_TOOLS),
           tool_b = factor(tool_b, levels = rev(CNV_TOOLS)))]
  ggplot(d, aes(tool_a, tool_b, fill = r)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 5 / .pt, family = FONT) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = .5,
                         limits = c(0, 1), breaks = c(0, .5, 1), name = "Pearson r") +
    coord_fixed() + labs(title = TOOL_TITLES[[key]], x = NULL, y = NULL) +
    guides(fill = if (key == "j") guide_colourbar() else "none") + theme_f5() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

## ==== geometry (inches, y from top) ===========================================
# Row 2 pairs the GT-correlation bar with the subclone correspondence; rows 3-4 pair
# one CNV heatmap with the matching spatial map per platform; row 5 is the 3x3 tool
# concordance, one square per platform.
FIG_W <- 7.09; Lm <- 0.06; G <- 0.12; RG <- 0.16; BG <- 0.24; TM <- 0.10; BM <- 0.06
xR <- FIG_W - Lm; Ha <- 1.13; Hc <- 1.95; Hp <- 1.82; Ht <- 1.55
r1 <- TM; r2 <- r1 + Ha + RG; r3 <- r2 + Hc + BG; r4 <- r3 + Hp + RG; r5 <- r4 + Hp + BG
FIG_H <- r5 + Ht + BM
colW <- (xR - Lm - G) / 2; xg <- Lm + colW + G
bW   <- 1.95                                   # correlation bar, narrow left of row 2
tW   <- (xR - Lm - 2 * G) / 3                  # three equal squares on row 5
B <- rbind(
  data.frame(k = "a", lab = "a", x = Lm, yt = r1, w = xR - Lm, h = Ha),
  data.frame(k = "b", lab = "b", x = Lm, yt = r2, w = bW, h = Hc),
  data.frame(k = "c", lab = "c", x = Lm + bW + G, yt = r2, w = xR - (Lm + bW + G), h = Hc),
  data.frame(k = "d", lab = "d", x = Lm, yt = r3, w = colW, h = Hp),
  data.frame(k = "e", lab = "e", x = xg, yt = r3, w = colW, h = Hp),
  data.frame(k = "f", lab = "f", x = Lm, yt = r4, w = colW, h = Hp),
  data.frame(k = "g", lab = "g", x = xg, yt = r4, w = colW, h = Hp),
  data.frame(k = "h", lab = "h", x = Lm, yt = r5, w = tW, h = Ht),
  data.frame(k = "i", lab = "i", x = Lm + tW + G, yt = r5, w = tW, h = Ht),
  data.frame(k = "j", lab = "j", x = Lm + 2 * (tW + G), yt = r5, w = tW, h = Ht))

message(sprintf("building fig5 ... page %.2f x %.2f in (%.0f x %.0f mm)", FIG_W, FIG_H, FIG_W * 25.4, FIG_H * 25.4))
asp <- function(k) { b <- B[B$k == k, ]; b$w / b$h }
P <- list(a = build_a(), b = build_b(), c = build_c(),
          d = build_d(), e = build_e(asp("e")), f = build_f(), g = build_g(asp("g")),
          h = build_tool_r("h"), i = build_tool_r("i"), j = build_tool_r("j"))

## ---- standalone per-panel exports ----
save_panel(P$a, "fig5_panel_a_schematic",               7.0, 1.7)
save_panel(P$b, "fig5_panel_b_pseudobulk_corr",         2.0, 1.9)
save_panel(P$c, "fig5_panel_c_subclone_correspondence", 5.0, 2.2)
save_panel(P$d, "fig5_panel_d_atera_subclone_heatmap",  3.4, 2.0)
save_panel(P$e, "fig5_panel_e_atera_subclone_spatial",  3.4, 2.0)
save_panel(P$f, "fig5_panel_f_strataMap_cnv_heatmap",     3.4, 2.0)
save_panel(P$g, "fig5_panel_g_strataMap_cnv_spatial",     3.4, 2.0)
save_panel(P$h, "fig5_panel_h_tool_concordance_visiumHD", 2.3, 1.7)
save_panel(P$i, "fig5_panel_i_tool_concordance_xenium",   2.3, 1.7)
save_panel(P$j, "fig5_panel_j_tool_concordance_strataMap",  2.9, 1.7)

## ---- compose ----
canvas <- ggdraw()
for (i in seq_len(nrow(B))) { b <- B[i, ]
  canvas <- canvas + draw_plot(P[[b$k]], x = b$x / FIG_W, y = (FIG_H - b$yt - b$h) / FIG_H,
                               width = b$w / FIG_W, height = b$h / FIG_H) }
canvas <- canvas + draw_plot_label(B$lab, x = B$x / FIG_W, y = (FIG_H - B$yt) / FIG_H,
                                   hjust = 0, vjust = 1, size = S_LET, family = FONT, fontface = "bold")
save_panel_both(canvas, "fig5_compose", OUTF, FIG_W, FIG_H, FONT)
message(sprintf("wrote fig5_compose.{pdf,png}  (%.0f x %.0f mm)", FIG_W * 25.4, FIG_H * 25.4))
