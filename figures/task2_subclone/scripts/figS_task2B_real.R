#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig (Task 2, REAL spatial-coherence diagnostics + examples) —
# rebuilt from scratch using ONLY task2_subclone/input_data/real.
#
# Panels (letters = layout order):
#   a predicted # subclones SAMPLE x tool heatmap        [bestcase.csv]
#   b assortativity SAMPLE x tool heatmap (mono/blank)   [bestcase.csv]
#   c pairwise paired-Wilcoxon on assortativity (FFPE|FF-BAM)  [bestcase.csv]
#   d assortativity z-score box FFPE|FF-BAM              [bestcase.csv]
#   e predicted subclones vs z-score                    [bestcase.csv]
#   f spatial example metastatic_prostate_ffpe / pt15_P2 [real/spatial/...]
#   g spatial example metastatic_prostate_ffpe / pt10_A3 [real/spatial/...]
#
# Rule 6: the pairwise significance heatmap (c) colours AND stars use the
# BH-adjusted p (p_adj_BH), not the raw p. All plotted values computed here.
# Sources ONLY the two utils scripts. Panels -> supp_figures/, data -> data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
  library(png); library(dplyr); library(stringr); library(scales)
})
if (requireNamespace("ggrastr", quietly = TRUE)) library(ggrastr) else
  rasterise <- function(x, ...) x

.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("figS_task2B_real.R")
})()
SCRIPTS_DIR <- dirname(.this_file)
TASK_DIR    <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS_DIR   <- normalizePath(file.path(TASK_DIR, "..", "utils"))
IN_REAL     <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG     <- file.path(TASK_DIR, "supp_figures")
OUT_TAB     <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS_DIR, "figure_setup.R"))
source(file.path(UTILS_DIR, "task2_subclone_setup.R"))
write_tab <- function(df, name) readr::write_csv(as.data.frame(df), file.path(OUT_TAB, name))

FONT <- pick_font()
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

BASE    <- 7
sz      <- function(x) x * 0.72
pt      <- function(x) x * 0.6
PT_SIZE <- 0.42
theme_c <- function() {
  theme_task2_paper(base_size = BASE) +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0.5),
          plot.subtitle = element_text(size = BASE - 1.5, color = "grey35"),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

REAL_CSV <- file.path(IN_REAL, "subclone_spatial_assortativity_summary.csv")
SPATIAL  <- file.path(IN_REAL, "spatial")
REF_FREE_TOOLS <- c("STARCH", "CalicoST", "TUSCAN")

t2             <- t2_load_real_bestcase(REAL_CSV)
t2_best        <- t2$best
t2_nc          <- t2$nc
t2_best_ab     <- t2$best_ab
t2_tool_levels <- t2$tool_levels
t2_tool_pal    <- t2$tool_pal

# shared y-order for the two SAMPLE x tool heatmaps (expression tools by
# descending median assortativity, then allele-only CalicoST/Numbat at bottom).
.t2_assort_rank <- t2_best[, .(m = median(assortativity_obs, na.rm = TRUE)),
                           by = tool][order(-m), as.character(tool)]
t2_heat_tool_order <- c(setdiff(.t2_assort_rank, c("CalicoST", "Numbat")),
                        intersect(.t2_assort_rank, c("CalicoST", "Numbat")))

# =====================================================================
# (d) z-score boxplot, faceted FFPE | FF-BAM
# =====================================================================
ord_d <- t2_facet_median_order(t2_best_ab, "z_score")
lab_d <- t2_group_omni(t2_best_ab, "z_score")
p_zbox <- ggplot(ord_d$d, aes(x = xid, y = z_score, fill = tool)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  geom_jitter(aes(color = tool), width = 0.15, size = pt(2.6), alpha = 0.8) +
  geom_text(data = lab_d, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.5, size = sz(3.0), fontface = "italic") +
  facet_wrap(~ group, scales = "free_x") +
  scale_x_discrete(labels = ord_d$labels) +
  scale_fill_manual(values = t2_tool_pal) + scale_color_manual(values = t2_tool_pal) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "Spatial assortativity z-score (FFPE | FF-BAM)",
       x = NULL, y = "Z-score (vs permutation null)") +
  theme_c() + theme(legend.position = "none")
write_tab(ord_d$d %>% as.data.frame() %>% transmute(group, tool, dataset, sample, z_score),
          "figS_task2B_panel_d_zscore.csv")
write_tab(lab_d, "figS_task2B_panel_d_omnibus.csv")

# =====================================================================
# (e) subclones vs z-score, faceted FFPE | FF-BAM (Spearman rho + p)
# =====================================================================
lab_e <- t2_best[is.finite(n_subclones) & is.finite(z_score), {
  ct <- suppressWarnings(cor.test(n_subclones, z_score, method = "spearman"))
  pl <- if (ct$p.value < 0.001) "italic(p) < 0.001" else sprintf("italic(p) == %.3f", ct$p.value)
  .(rho = unname(ct$estimate), p_value = ct$p.value,
    label = sprintf("italic(rho)[s] == %.2f * \",\" ~ %s", ct$estimate, pl))
}, by = group]
subclone_x_log <- list(
  scale_x_log10(breaks = c(2, 3, 5, 10, 20, 50), minor_breaks = NULL,
                expand = expansion(mult = c(0.04, 0.04))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(1, "pt"), mid = unit(2, "pt"), long = unit(3, "pt")),
  labs(x = "Predicted # subclones (log10)"))
p_svz <- ggplot(t2_best, aes(x = n_subclones, y = z_score, color = tool)) +
  geom_point(size = pt(2.4), alpha = 0.8) +
  geom_text(data = lab_e, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.6, size = sz(3.0), parse = TRUE) +
  facet_wrap(~ group) +
  scale_color_manual(values = t2_tool_pal, name = "Tool") +
  subclone_x_log +
  labs(title = "Subclones vs z-score", y = "Z-score") +
  theme_c() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.size = unit(0.22, "cm"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2))
write_tab(t2_best %>% as.data.frame() %>% transmute(group, tool, dataset, sample, n_subclones, z_score),
          "figS_task2B_panel_e_subclones_vs_zscore.csv")
write_tab(lab_e %>% as.data.frame() %>% transmute(group, spearman_rho = rho, p_value),
          "figS_task2B_panel_e_spearman.csv")

# =====================================================================
# (b) assortativity heatmap SAMPLE x tool (monoclone / blank)
# =====================================================================
dnc <- copy(t2_nc)
dnc[, samp_key := paste(dataset, sample, sep = "::")]
dnc[, cell_val := fifelse(status == "ok" & is.finite(assortativity_obs), assortativity_obs, NA_real_)]
dnc[, is_mono  := status == "single_subclone" | (is.finite(n_subclones) & n_subclones == 1)]
d_key_order <- unique(dnc[order(group, dataset, sample)]$samp_key)
samp_dt <- unique(dnc[, .(samp_key, sample, group)])
grid <- samp_dt[, .(tool = t2_tool_levels), by = .(samp_key, sample, group)]
hd <- merge(grid, dnc[, .(samp_key, tool = as.character(tool), cell_val, is_mono)],
            by = c("samp_key", "tool"), all.x = TRUE)
hd[is.na(is_mono), is_mono := FALSE]
hd[, cell_lab := fifelse(!is.na(cell_val), sprintf("%.2f", cell_val),
                 fifelse(is_mono, "monoclone", ""))]
hd[, lab_size := fifelse(cell_lab == "monoclone", 1.2, 1.8)]
hd[, samp_key := factor(samp_key, levels = d_key_order)]
hd[, tool := factor(as.character(tool), levels = rev(t2_heat_tool_order))]
d_samp_lab <- setNames(samp_dt$sample, samp_dt$samp_key)
p_assort_hm <- ggplot(hd, aes(x = samp_key, y = tool, fill = cell_val)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_lab, size = lab_size)) +
  scale_size_identity() +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = d_samp_lab) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
                       midpoint = 0, na.value = "grey90", name = "Assortativity") +
  labs(title = "Assortativity by sample × tool ('monoclone' = single-clone; blank = not run)",
       x = NULL, y = NULL) +
  theme_c() +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.30, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 3),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2))
write_tab(hd %>% as.data.frame() %>%
            transmute(samp_key, sample = d_samp_lab[as.character(samp_key)], tool,
                      assortativity = cell_val, is_mono, cell_lab),
          "figS_task2B_panel_b_assortativity_heatmap.csv")

# =====================================================================
# (a) predicted # subclones heatmap SAMPLE x tool
# =====================================================================
ncc <- copy(t2_nc)
ncc[, samp_key := paste(dataset, sample, sep = "::")]
key_order <- unique(ncc[order(group, dataset, sample)]$samp_key)
wnc <- dcast(ncc, samp_key + sample + group ~ tool, value.var = "n_subclones")
hn  <- melt(wnc, id.vars = c("samp_key", "sample", "group"),
            variable.name = "tool", value.name = "n_subclones")
hn[, samp_key := factor(samp_key, levels = key_order)]
hn[, tool := factor(as.character(tool), levels = rev(t2_heat_tool_order))]
n_samp_lab <- setNames(unique(ncc[, .(samp_key, sample)])$sample,
                       unique(ncc[, .(samp_key, sample)])$samp_key)
p_ncount_hm <- ggplot(hn, aes(x = samp_key, y = tool, fill = n_subclones)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = fifelse(is.na(n_subclones), "", as.character(n_subclones))), size = 1.8) +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = n_samp_lab) +
  scale_fill_gradient(low = "#FFF7BC", high = "#800026", trans = "log10",
                      breaks = c(1, 2, 5, 10, 20, 50), na.value = "grey90", name = "# subclones") +
  labs(title = "Predicted # subclones per sample × tool (blank = not run)", x = NULL, y = NULL) +
  theme_c() +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.30, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 3),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2))
write_tab(hn %>% as.data.frame() %>%
            transmute(samp_key, sample = n_samp_lab[as.character(samp_key)], group, tool, n_subclones),
          "figS_task2B_panel_a_nsubclones_heatmap.csv")

# =====================================================================
# (c) pairwise paired-Wilcoxon on assortativity, per group.
#     Rule 6: fill = -log10(p_adj_BH); stars from p_adj_BH.
# =====================================================================
pw_assort <- t2_group_pairwise(t2_best, "assortativity_obs")
write_tab(pw_assort %>% as.data.frame() %>%
            transmute(platform, group1, group2, n_pairs, median1, median2,
                      statistic, p_value, p_adj_BH),
          "figS_task2B_panel_c_pairwise_wilcox.csv")
sig_star <- function(p) fifelse(is.na(p), "",
            fifelse(p <= 0.001, "***", fifelse(p <= 0.01, "**", fifelse(p <= 0.05, "*", ""))))
make_sym <- function(sub, tools) {
  d1 <- sub[, .(row = group1, col = group2, p_adj_BH)]
  d2 <- sub[, .(row = group2, col = group1, p_adj_BH)]
  dg <- data.table(row = tools, col = tools, p_adj_BH = NA_real_)
  h  <- rbindlist(list(d1, d2, dg), use.names = TRUE)
  h[, star := sig_star(p_adj_BH)]
  h[, nlp  := fifelse(is.finite(p_adj_BH), -log10(pmax(p_adj_BH, .Machine$double.xmin)), NA_real_)]
  h[, row  := factor(row, levels = rev(tools))]
  h[, col  := factor(col, levels = tools)]
  h
}
nlp_max <- max(-log10(pmax(pw_assort$p_adj_BH, .Machine$double.xmin)), na.rm = TRUE)
pw_fill_limits <- c(0, nlp_max)
make_pw_plot <- function(grp, cap = NULL) {
  sub   <- pw_assort[platform == grp]
  tools <- intersect(t2_tool_levels, unique(c(sub$group1, sub$group2)))
  ggplot(make_sym(sub, tools), aes(col, row, fill = nlp)) +
    geom_tile(color = "grey88", linewidth = 0.22) +
    geom_text(aes(label = star), size = 2.0, vjust = 0.75, color = "grey15", na.rm = TRUE) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = expression(-log[10](p[adj])),
                        na.value = "grey95", limits = pw_fill_limits, oob = scales::squish) +
    coord_equal() +
    labs(title = paste0("Pairwise Wilcoxon · ", grp), x = NULL, y = NULL, caption = cap) +
    theme_pub(base_size = BASE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 2),
          axis.text.y = element_text(size = BASE - 2),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.25, "cm"), plot.caption = element_text(size = BASE - 2, hjust = 1),
          plot.title = element_text(size = BASE, face = "bold", hjust = 0.5))
}
p_pw <- ((make_pw_plot("FFPE") + labs(tag = "c")) |
         make_pw_plot("FF-BAM", cap = "BH: * .05  ** .01  *** .001")) +
  plot_layout(guides = "collect")

# =====================================================================
# SPATIAL panels f / g (read the clean input_data/real/spatial layout)
# =====================================================================
norm_bc <- function(x) gsub("\\.", "-", trimws(as.character(x)))
load_positions <- function(dir) {
  pos <- fread(file.path(dir, "positions.csv")); setnames(pos, 1, "barcode")
  pos[, .(barcode = norm_bc(barcode),
          x = as.numeric(pxl_col_in_fullres), y = as.numeric(pxl_row_in_fullres))]
}
load_mask <- function(dir) {
  rd <- function(f) { p <- file.path(dir, f)
    if (!file.exists(p)) return(character(0)); b <- norm_bc(readLines(p, warn = FALSE)); b[nzchar(b)] }
  setdiff(unique(rd("keep_barcodes.txt")), unique(rd("ref_barcodes.txt")))
}
he_subpanel <- function(dir, he_title, tag = NULL) {
  img <- readPNG(file.path(dir, "he_lowres.png")); h <- dim(img)[1]; w <- dim(img)[2]
  right_pad <- if (h / w > 1) (h / w - 1) else 0
  p <- ggplot() +
    annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h) +
    scale_x_continuous(limits = c(0, w), expand = expansion(mult = c(0, right_pad))) +
    scale_y_continuous(limits = c(0, h), expand = c(0, 0), trans = "reverse") +
    labs(title = paste0("H&E\n", he_title)) +
    theme_void(base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = BASE, hjust = 0, lineheight = 0.9),
          plot.margin = margin(1, 1, 1, 1))
  if (!is.null(tag)) p <- p + labs(tag = tag); p
}
tool_subpanel <- function(dir, tool, pos, mask, ds, samp) {
  pf <- file.path(dir, paste0("subclone__", tool, ".csv"))
  if (!file.exists(pf))
    return(list(plot = ggplot() + labs(title = sprintf("%s\n(missing)", tool_label(tool))) +
                  theme_void(base_family = FONT) +
                  theme(plot.title = element_text(face = "bold", size = BASE, hjust = 0.5, lineheight = 0.9)),
                pts = NULL, summ = data.table(sample = samp, tool = tool, k_clones = NA_integer_,
                                              assortativity_r = NA_real_, n_mask_spots = length(mask))))
  pred <- fread(pf)
  pred <- pred[, .(barcode = norm_bc(spot_id), subclone = as.character(subclone))]
  pred <- unique(pred[!is.na(barcode) & nzchar(barcode) & !is.na(subclone)], by = "barcode")
  d <- merge(pos, pred, by = "barcode", all.x = TRUE)
  d[, in_tumor := barcode %in% mask]
  non_tumor  <- d[in_tumor == FALSE]
  unassigned <- d[in_tumor == TRUE & is.na(subclone)]
  fg         <- d[in_tumor == TRUE & !is.na(subclone)]
  lev <- fg[, .N, by = subclone][order(-N), subclone]
  fg[, subclone := factor(subclone, levels = lev)]; setorder(fg, subclone)
  k <- length(lev)
  pal <- setNames(scales::hue_pal()(max(k, 1))[seq_len(k)], lev)
  xr <- range(pos$x, na.rm = TRUE); yr <- range(pos$y, na.rm = TRUE)
  cx <- mean(xr); cy <- mean(yr); half <- max(diff(xr), diff(yr)) / 2 * 1.04
  r <- t2$assort(ds, samp, tool)
  title <- if (is.finite(r)) sprintf("%s\n%d clones · r=%.2f", tool_label(tool), k, r)
           else sprintf("%s\n(%d clones)", tool_label(tool), k)
  p <- ggplot() +
    rasterise(geom_point(data = non_tumor, aes(x, y), color = "grey88", size = PT_SIZE, alpha = 0.7, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = unassigned, aes(x, y), color = "grey40", size = PT_SIZE, alpha = 0.9, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = fg, aes(x, y, color = subclone), size = PT_SIZE, alpha = 0.95, stroke = 0), dpi = 300) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(cx - half, cx + half), expand = c(0, 0)) +
    scale_y_continuous(limits = c(cy - half, cy + half), expand = c(0, 0), trans = "reverse") +
    labs(title = title) +
    theme_void(base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = BASE, hjust = 0.5, lineheight = 0.9),
          plot.margin = margin(1, 1, 1, 1))
  pts <- rbind(
    non_tumor[,  .(sample = samp, tool = tool, barcode, x, y, category = "non_tumor",  subclone = NA_character_)],
    unassigned[, .(sample = samp, tool = tool, barcode, x, y, category = "unassigned", subclone = NA_character_)],
    fg[,         .(sample = samp, tool = tool, barcode, x, y, category = "subclone",   subclone = as.character(subclone))])
  list(plot = p, pts = pts,
       summ = data.table(sample = samp, tool = tool, k_clones = k,
                         assortativity_r = r, n_mask_spots = length(mask)))
}
grey_legend_h <- function() {
  d <- data.frame(x = c(1, 2.2), y = c(1, 1), lab = c("non-tumor", "unassigned"),
                  col = c("grey88", "grey40"))
  ggplot(d, aes(x, y)) +
    geom_point(aes(fill = col), shape = 21, color = "grey55", size = 1.6, stroke = 0.3) +
    geom_text(aes(label = lab), hjust = 0, nudge_x = 0.10, size = sz(2.6)) +
    scale_fill_identity() + scale_x_continuous(limits = c(0.8, 4.0)) +
    scale_y_continuous(limits = c(0.5, 1.5)) + theme_void(base_family = FONT)
}

TOOLS_SPATIAL <- c("InferCNV", "SCEVAN", "STARCH", "TUSCAN")
build_spatial_row <- function(ds, samp, he_title, tag) {
  dir  <- file.path(SPATIAL, ds, samp)
  pos  <- load_positions(dir); mask <- load_mask(dir)
  he   <- he_subpanel(dir, he_title, tag = tag)
  subs <- lapply(TOOLS_SPATIAL, function(t) tool_subpanel(dir, t, pos, mask, ds, samp))
  list(row = c(list(he), lapply(subs, `[[`, "plot")), subs = subs)
}
row_f <- build_spatial_row("metastatic_prostate_ffpe", "pt15_P2", "pt15_P2 (FFPE)", "f")
row_g <- build_spatial_row("metastatic_prostate_ffpe", "pt10_A3", "pt10_A3 (FFPE)", "g")
ef_block <- wrap_plots(c(row_f$row, row_g$row), ncol = 5, nrow = 2)

spatial_pts  <- rbindlist(c(lapply(row_f$subs, `[[`, "pts"), lapply(row_g$subs, `[[`, "pts")), fill = TRUE)
spatial_summ <- rbindlist(c(lapply(row_f$subs, `[[`, "summ"), lapply(row_g$subs, `[[`, "summ")), fill = TRUE)
write_tab(spatial_pts[sample == "pt15_P2"], "figS_task2B_panel_f_spatial_points.csv")
write_tab(spatial_summ[sample == "pt15_P2"], "figS_task2B_panel_f_spatial_summary.csv")
write_tab(spatial_pts[sample == "pt10_A3"], "figS_task2B_panel_g_spatial_points.csv")
write_tab(spatial_summ[sample == "pt10_A3"], "figS_task2B_panel_g_spatial_summary.csv")

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
row_de <- (p_zbox + labs(tag = "d")) | (p_svz + labs(tag = "e"))
final <- (p_ncount_hm + labs(tag = "a")) / (p_assort_hm + labs(tag = "b")) / p_pw / row_de /
  ef_block / grey_legend_h() +
  plot_layout(heights = c(1.5, 1.5, 2.4, 2.3, 3.1, 0.28)) &
  theme(plot.tag = element_text(face = "bold", size = 9), text = element_text(family = FONT))
W <- 7.09; H <- 11.8
ggsave(file.path(OUT_FIG, "figS_task2B_real_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "figS_task2B_real_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")

panel_f <- wrap_plots(row_f$row, nrow = 1, ncol = 5)
panel_g <- wrap_plots(row_g$row, nrow = 1, ncol = 5)
panels <- list(a = list(p_ncount_hm, 7.09, 2.0), b = list(p_assort_hm, 7.09, 2.0),
               c = list(p_pw, 5.0, 2.6), d = list(p_zbox, 4.2, 2.8), e = list(p_svz, 4.6, 2.8),
               f = list(panel_f, 7.09, 1.7), g = list(panel_g, 7.09, 1.7))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("figS_task2B_panel_", nm))
  pp <- pl[[1]] & theme(text = element_text(family = FONT))
  ggsave(paste0(base, ".pdf"), pp, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pp, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] figS_task2B_real_compose + panels a-g (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
cat("[Info] d z-score omnibus:\n"); print(lab_d)
cat("[Info] e Spearman (subclones vs z):\n"); print(lab_e[, .(group, rho, p_value)])
cat("[Info] c pairwise (BH) assortativity:\n"); print(pw_assort[, .(platform, group1, group2, n_pairs, p_value, p_adj_BH)])
cat("[Info] spatial summary:\n"); print(spatial_summ)
