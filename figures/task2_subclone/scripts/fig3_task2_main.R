#!/usr/bin/env Rscript
# =====================================================================
# Fig 3 (Task 2: subclone spatial concordance) — MAIN composed figure,
# rebuilt from scratch using ONLY task2_subclone/input_data.
#
# Panels (letters = reading order):
#   a multiclonal ARI boxplot (GT>=2)              [synth/subclone_ari_summary.csv]
#   b monoclonal specificity dotplot (truth=1)     [synth/subclone_ari_summary.csv]
#   c predicted vs true subclone count             [synth/subclone_ari_summary.csv]
#   d real spatial assortativity box FFPE|FF-BAM   [real/...bestcase.csv]
#   e predicted subclones vs assortativity         [real/...bestcase.csv]
#   f spatial example prostate_eva_ffpe / P10_P3   [real/spatial/...]
#   g spatial example nature_prostate_ff / H1_4    [real/spatial/...]
#
# Every panel is replotted from data (no PNG copy). All plotted values are
# computed here from input_data. Sources ONLY the two utils scripts. Emits each
# panel as pdf+png in main_figures/ and each panel's plotted data in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
  library(png); library(dplyr); library(stringr); library(scales)
})
if (requireNamespace("ggrastr", quietly = TRUE)) library(ggrastr) else
  rasterise <- function(x, ...) x

# ---- locate self / paths --------------------------------------------
.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("fig3_task2_main.R")
})()
SCRIPTS_DIR <- dirname(.this_file)
TASK_DIR    <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS_DIR   <- normalizePath(file.path(TASK_DIR, "..", "utils"))
IN_SYN      <- file.path(TASK_DIR, "input_data", "synth")
IN_REAL     <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG     <- file.path(TASK_DIR, "main_figures")
OUT_TAB     <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS_DIR, "figure_setup.R"))
source(file.path(UTILS_DIR, "task2_subclone_setup.R"))
write_tab <- function(df, name) readr::write_csv(as.data.frame(df), file.path(OUT_TAB, name))

FONT <- pick_font()
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

SYN_CSV  <- file.path(IN_SYN,  "subclone_ari_summary.csv")
REAL_CSV <- file.path(IN_REAL, "subclone_spatial_assortativity_summary.csv")
SPATIAL  <- file.path(IN_REAL, "spatial")

# ---- shared sizing (from the frozen fig3_compose.R) -----------------
BASE    <- 7
sz      <- function(x) x * 0.62
pt      <- function(x) x * 0.55
PT_SIZE <- 0.3

theme_c <- function() {
  theme_task2_paper(base_size = BASE) +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0.5),
          plot.caption = element_text(size = BASE - 2, color = "gray40", hjust = 1),
          axis.text.x = element_text(angle = 30, hjust = 1))
}
pval_tr <- function(label, size)
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.05, vjust = 1.5,
           size = size, fontface = "italic")
# display-only: append "*" to CalicoST / STARCH labels (defined in legend).
star_lab <- function(x) {
  v <- ifelse(as.character(x) %in% c("CalicoST", "STARCH"),
              paste0(as.character(x), "*"), as.character(x))
  if (!is.null(names(x))) names(v) <- names(x)
  v
}

# =====================================================================
# SYNTHETIC panels a / b / c
# =====================================================================
df <- load_subclone_summary(SYN_CSV)

# ---- (a) multiclonal ARI boxplot ------------------------------------
multi_ok  <- df %>% filter(status == "ok", n_gt_subclones >= 2)
multi_all <- df %>% filter(n_gt_subclones >= 2)
om_a <- block_omnibus(data.frame(value = multi_all$ari,
                                 group = multi_all$tool, block = multi_all$rep))
lab_a <- omnibus_label(om_a)
p_a <- ggplot(multi_ok, aes(x = tool, y = ari)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_boxplot(aes(fill = tool), alpha = 0.45, outlier.shape = NA,
               width = 0.55, color = "black", linewidth = 0.4) +
  geom_jitter(aes(fill = tool), width = 0.10, size = pt(2.4), alpha = 0.9,
              shape = 21, color = "white", stroke = 0.4) +
  pval_tr(lab_a, sz(3.6)) +
  scale_fill_manual(values = TOOL_COLORS_SYNTH_TASK2) +
  scale_x_discrete(labels = star_lab) +
  scale_y_continuous(limits = c(-0.20, 1.22), breaks = seq(0, 1, 0.2)) +
  labs(title = paste0("Multiclonal partitioning\naccuracy (n=",
                      length(unique(multi_ok$rep)), ")"),
       x = NULL, y = "Adjusted Rand Index (ARI)") +
  theme_c() + theme(legend.position = "none")
write_tab(multi_ok %>% transmute(tool, rep, n_gt_subclones, ari),
          "fig3_panel_a_multiclonal_ari.csv")
write_tab(omnibus_row(om_a, "multiclonal_ari"), "fig3_panel_a_omnibus.csv")

# ---- (b) monoclonal specificity dotplot -----------------------------
mono_df <- df %>% filter(n_gt_subclones == 1) %>%
  mutate(is_failure  = (status == "missing_prediction"),
         n_pred_plot = ifelse(is_failure, NA_real_, n_pred_subclones))
om_b <- block_omnibus(data.frame(value = mono_df$n_pred_plot,
                                 group = mono_df$tool, block = mono_df$rep))
lab_b  <- omnibus_label(om_b)
n_mono <- length(unique(mono_df$rep))
p_b <- ggplot(mono_df, aes(x = tool, y = n_pred_plot)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", linewidth = 0.5) +
  geom_jitter(aes(fill = tool), width = 0.12, size = pt(3.8),
              shape = 21, color = "white", stroke = 0.5, alpha = 0.92, na.rm = TRUE) +
  geom_point(data = filter(mono_df, is_failure, tool != "Numbat"), aes(x = tool, y = 0.5),
             shape = 4, size = pt(4), color = "gray45", stroke = 0.9, inherit.aes = FALSE) +
  annotate("text", x = 0.55, y = 1.35, label = "Truth = 1",
           color = "darkgreen", size = sz(3), hjust = 0, fontface = "italic") +
  pval_tr(lab_b, sz(3.6)) +
  scale_fill_manual(values = TOOL_COLORS_SYNTH_TASK2) +
  scale_x_discrete(labels = star_lab) +
  scale_y_continuous(trans = "log1p", breaks = c(1, 2, 5, 10, 20),
                     labels = c("1", "2", "5", "10", "20"), limits = c(0.3, 25)) +
  labs(title = paste0("Monoclonal specificity\n(truth = 1, n=", n_mono, ")"),
       x = NULL, y = "Predicted # subclones (log)") +
  theme_c() + theme(legend.position = "none")
write_tab(mono_df %>% transmute(tool, rep, n_gt_subclones, n_pred_subclones, status, is_failure),
          "fig3_panel_b_monoclonal_specificity.csv")
write_tab(omnibus_row(om_b, "monoclonal_specificity"), "fig3_panel_b_omnibus.csv")

# ---- (c) predicted vs true subclone count ---------------------------
ok_df <- df %>% filter(status == "ok")
p_c <- ggplot(ok_df, aes(x = n_gt_subclones, y = n_pred_subclones)) +
  annotate("rect", xmin = 0.3, xmax = 4.7, ymin = 8, ymax = 26, alpha = 0.08, fill = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_jitter(aes(fill = tool), width = 0.12, height = 0.12, size = pt(2.8), alpha = 0.78,
              shape = 21, color = "white", stroke = 0.4) +
  annotate("text", x = 4.55, y = 23, label = "Severe over-\nsegmentation",
           color = "#B53C2A", fontface = "bold", size = sz(3), hjust = 1) +
  annotate("text", x = 1, y = 1, label = "y = x",
           color = "gray40", size = sz(3), hjust = -0.3, vjust = 1.5) +
  scale_fill_manual(values = TOOL_COLORS_SYNTH_TASK2, name = "Tool", labels = star_lab) +
  scale_x_continuous(breaks = 1:4, limits = c(0.3, 4.7)) +
  scale_y_continuous(limits = c(-0.5, 26)) +
  labs(title = "Predicted vs true\nsubclone count",
       x = "Ground truth # subclones", y = "Predicted # subclones") +
  theme_c() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.size = unit(0.22, "cm"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.3))
write_tab(ok_df %>% transmute(tool, rep, n_gt_subclones, n_pred_subclones),
          "fig3_panel_c_pred_vs_true_count.csv")

# =====================================================================
# REAL panels d + e (FFPE | FF-BAM), facets stacked vertically
# =====================================================================
REF_FREE_TOOLS <- c("STARCH", "CalicoST", "TUSCAN")
t2          <- t2_load_real_bestcase(REAL_CSV)
t2_best     <- t2$best
t2_best_ab  <- t2$best_ab
t2_tool_pal <- t2$tool_pal

# ---- (d) assortativity boxplot faceted FFPE | FF-BAM ----------------
ord_d <- t2_facet_median_order(t2_best_ab, "assortativity_obs")
lab_d <- t2_group_omni(t2_best_ab, "assortativity_obs")
p_d <- ggplot(ord_d$d, aes(x = xid, y = assortativity_obs, fill = tool)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  geom_jitter(aes(color = tool), width = 0.15, size = pt(2.6), alpha = 0.8) +
  geom_text(data = lab_d, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.5, size = sz(3.4), fontface = "italic") +
  facet_wrap(~ group, scales = "free_x", ncol = 1) +
  scale_x_discrete(labels = star_lab(ord_d$labels)) +
  scale_fill_manual(values = t2_tool_pal) + scale_color_manual(values = t2_tool_pal) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "Real-data spatial\nassortativity", x = NULL, y = "Assortativity") +
  theme_c() + theme(legend.position = "none")
write_tab(ord_d$d %>% as.data.frame() %>%
            transmute(group, tool, dataset, sample, assortativity_obs),
          "fig3_panel_d_assortativity.csv")
write_tab(lab_d, "fig3_panel_d_omnibus.csv")

# ---- (e) predicted subclones vs assortativity, faceted --------------
subclone_x_log <- list(
  scale_x_log10(breaks = c(2, 3, 5, 10, 20, 50), minor_breaks = NULL,
                expand = expansion(mult = c(0.04, 0.04))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(1, "pt"), mid = unit(2, "pt"), long = unit(3, "pt")),
  labs(x = "Predicted # subclones (log10)"))
cor_e <- t2_best[, {
  ct   <- suppressWarnings(cor.test(n_subclones, assortativity_obs, method = "spearman"))
  pstr <- scales::pvalue(ct$p.value, accuracy = 0.001)
  psym <- if (startsWith(pstr, "<")) sprintf('italic(p) < "%s"', sub("^<", "", pstr))
          else                       sprintf('italic(p) == "%s"', pstr)
  .(rho = unname(ct$estimate), p_value = ct$p.value,
    label = sprintf("rho == %.2f * ',' ~ %s", unname(ct$estimate), psym))
}, by = group]
p_e <- ggplot(t2_best, aes(x = n_subclones, y = assortativity_obs, color = tool)) +
  geom_point(size = pt(2.4), alpha = 0.8) +
  geom_text(data = cor_e, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.5, size = sz(3.4), parse = TRUE) +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = t2_tool_pal, name = "Tool", labels = star_lab) +
  subclone_x_log +
  labs(title = "Subclones vs\nassortativity", y = "Assortativity") +
  theme_c() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.size = unit(0.22, "cm"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2))
write_tab(t2_best %>% as.data.frame() %>%
            transmute(group, tool, dataset, sample, n_subclones, assortativity_obs),
          "fig3_panel_e_subclones_vs_assortativity.csv")
write_tab(cor_e %>% as.data.frame() %>% transmute(group, spearman_rho = rho, p_value),
          "fig3_panel_e_spearman.csv")

# =====================================================================
# SPATIAL panels f / g  (read the clean input_data/real/spatial layout)
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
# Returns list(plot, table rows, summary row) for one sample x tool scatter.
tool_subpanel <- function(dir, tool, pos, mask, ds, samp) {
  pf <- file.path(dir, paste0("subclone__", tool, ".csv"))
  if (!file.exists(pf))
    return(list(plot = ggplot() + labs(title = sprintf("%s\n(missing)", star_lab(tool_label(tool)))) +
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
  title <- if (is.finite(r)) sprintf("%s\n%d clones · r=%.2f", star_lab(tool_label(tool)), k, r)
           else sprintf("%s\n(%d clones)", star_lab(tool_label(tool)), k)
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
grey_legend_v <- function() {
  d <- data.frame(x = c(1, 1), y = c(2, 1), lab = c("non-tumor", "unassigned"),
                  col = c("grey88", "grey40"))
  ggplot(d, aes(x, y)) +
    geom_point(aes(fill = col), shape = 21, color = "grey55", size = 1.6, stroke = 0.3) +
    geom_text(aes(label = lab), hjust = 0, nudge_x = 0.12, size = sz(2.8)) +
    scale_fill_identity() + scale_x_continuous(limits = c(0.8, 3.0)) +
    scale_y_continuous(limits = c(-1, 4)) + theme_void(base_family = FONT)
}

# ---- panel f: prostate_eva_ffpe / P10_P3 ----------------------------
f_dir  <- file.path(SPATIAL, "prostate_eva_ffpe", "P10_P3")
f_pos  <- load_positions(f_dir); f_mask <- load_mask(f_dir)
f_tools <- c("InferCNV", "SCEVAN", "STARCH", "TUSCAN")
he_f   <- he_subpanel(f_dir, "P10_P3 (FFPE)")
sub_f  <- lapply(f_tools, function(t) tool_subpanel(f_dir, t, f_pos, f_mask, "prostate_eva_ffpe", "P10_P3"))

# ---- panel g: nature_prostate_ff / H1_4 -----------------------------
g_dir  <- file.path(SPATIAL, "nature_prostate_ff", "H1_4")
g_pos  <- load_positions(g_dir); g_mask <- load_mask(g_dir)
g_expr <- c("InferCNV", "SCEVAN", "STARCH", "TUSCAN")
g_all  <- c("CalicoST", "Numbat")
he_g   <- he_subpanel(g_dir, "H1_4 (FF)")
sub_g_expr <- lapply(g_expr, function(t) tool_subpanel(g_dir, t, g_pos, g_mask, "nature_prostate_ff", "H1_4"))
sub_g_all  <- lapply(g_all,  function(t) tool_subpanel(g_dir, t, g_pos, g_mask, "nature_prostate_ff", "H1_4"))

# export spatial data + summaries
spatial_pts  <- rbindlist(c(lapply(sub_f, `[[`, "pts"),
                            lapply(sub_g_expr, `[[`, "pts"), lapply(sub_g_all, `[[`, "pts")), fill = TRUE)
spatial_summ <- rbindlist(c(lapply(sub_f, `[[`, "summ"),
                            lapply(sub_g_expr, `[[`, "summ"), lapply(sub_g_all, `[[`, "summ")), fill = TRUE)
write_tab(spatial_pts[sample == "P10_P3"], "fig3_panel_f_spatial_points.csv")
write_tab(spatial_summ[sample == "P10_P3"], "fig3_panel_f_spatial_summary.csv")
write_tab(spatial_pts[sample == "H1_4"], "fig3_panel_g_spatial_points.csv")
write_tab(spatial_summ[sample == "H1_4"], "fig3_panel_g_spatial_summary.csv")

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
row_1  <- (p_a + labs(tag = "a")) | (p_b + labs(tag = "b")) | (p_c + labs(tag = "c"))
row_de <- (p_d + labs(tag = "d")) | (p_e + labs(tag = "e"))
ef_block <- wrap_plots(
  c(list(he_subpanel(f_dir, "P10_P3 (FFPE)", tag = "f")), lapply(sub_f, `[[`, "plot"),
    list(he_subpanel(g_dir, "H1_4 (FF)", tag = "g")), lapply(sub_g_expr, `[[`, "plot"),
    list(plot_spacer()), lapply(sub_g_all, `[[`, "plot"), list(grey_legend_v(), plot_spacer())),
  ncol = 5, nrow = 3)
final <- row_1 / row_de / ef_block +
  plot_layout(heights = c(1.8, 3.8, 5.25)) &
  theme(plot.tag = element_text(face = "bold", size = 9), text = element_text(family = FONT))

W <- 7.09; H <- 10.8
ggsave(file.path(OUT_FIG, "fig3_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig3_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")

# per-panel: a-e are single plots; f/g are the spatial rows for one sample.
panel_f <- wrap_plots(c(list(he_subpanel(f_dir, "P10_P3 (FFPE)")), lapply(sub_f, `[[`, "plot")),
                      nrow = 1, ncol = 5)
panel_g <- wrap_plots(c(list(he_subpanel(g_dir, "H1_4 (FF)")), lapply(sub_g_expr, `[[`, "plot"),
                        list(plot_spacer()), lapply(sub_g_all, `[[`, "plot"),
                        list(grey_legend_v(), plot_spacer())), ncol = 5, nrow = 2)
panels <- list(a = list(p_a, 2.5, 2.6), b = list(p_b, 2.5, 2.6), c = list(p_c, 3.2, 2.6),
               d = list(p_d, 3.4, 4.0), e = list(p_e, 3.6, 4.0),
               f = list(panel_f, 7.09, 1.9), g = list(panel_g, 7.09, 3.4))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("fig3_panel_", nm))
  pp <- pl[[1]] & theme(text = element_text(family = FONT))
  ggsave(paste0(base, ".pdf"), pp, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pp, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] fig3_compose + panels a-g (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
cat(sprintf("[Info] a omnibus: %s | b omnibus: %s (n_mono=%d)\n", lab_a, lab_b, n_mono))
cat("[Info] d assortativity omnibus:\n"); print(lab_d)
cat("[Info] e Spearman:\n"); print(cor_e[, .(group, rho, p_value)])
cat(sprintf("[Info] masks: P10_P3=%d  H1_4=%d spots\n", length(f_mask), length(g_mask)))
