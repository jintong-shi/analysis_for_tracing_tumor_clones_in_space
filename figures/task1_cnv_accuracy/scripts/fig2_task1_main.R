#!/usr/bin/env Rscript
# =====================================================================
# Fig 2 (Task 1: CNV-profile prediction) — MAIN figure, rebuilt from scratch.
#
# Reproduces paper/figures/fig2_merge/fig2_compose.R using ONLY the curated
# input_data of figures_final/task1_cnv_accuracy. Every plotted value is
# recomputed here; nothing is copied from a rendered figure.
#
# Panels:
#   a overall per-tool accuracy (Correlation | Max F1)         [synth normalref]
#   b broad vs narrow CNVs (Max F1)                            [synth normalref]
#   c per-tool CNV accuracy vs tumor purity (Max F1)           [synth purity]
#   d accuracy across datasets (Max F1)                        [real BEST]
#   e Ding FFPE per tool (Correlation, Max F1)                 [real BEST]
#   f Ding FF per tool (Correlation, Max F1)                   [real BEST]
#   g gene-level CNV heatmap, HT397B1-H2U1 (WES + tools)       [real WES + pseudobulk_gene]
#
# The ONLY script sourced is utils/figure_setup.R (shared settings + theme +
# Skillings-Mack). Outputs: composed fig + each panel (pdf/png) in main_figures/,
# per-panel plotted data in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(readr); library(purrr); library(ggrepel); library(tibble)
})

# ---- paths (this script's location -> project tree) -----------------
.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("fig2_task1_main.R")
})()
SCRIPTS_DIR <- dirname(.this_file)
TASK_DIR    <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS       <- normalizePath(file.path(TASK_DIR, "..", "utils", "figure_setup.R"))
IN_SYN      <- file.path(TASK_DIR, "input_data", "synth")
IN_REAL     <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG     <- file.path(TASK_DIR, "main_figures")
OUT_TAB     <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(UTILS)
source(file.path(.figure_setup_dir, "final_record_select.R"))

FONT <- pick_font()   # real "Arial" if installed, else bundled "Liberation Sans" (Arial metric clone), else "sans"
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

BASE_SIZE    <- 7
METRIC_ALL   <- c("spearman", "max_f1", "pauc_gain", "pauc_loss")
SCOPE_LEVELS <- c("overall", "broad_only", "narrow_only")
P_ANNO_SIZE  <- 5.5 / .pt
P_ANNO_FACE  <- "plain"

t1 <- figure_settings("tool_order", "synth_task1")
TOOL_ORDER_T1  <- t1$label
TOOL_COLORS_T1 <- setNames(t1$color, t1$label)
ID2LABEL_T1    <- setNames(t1$label, t1$id)
METRIC_LABELS  <- setting_vector("metric", "global", "label")
METRIC_LABELS["spearman"] <- "Correlation"

fmt_p_disp <- function(p) ifelse(is.na(p), "p = NA",
                                 ifelse(p < 0.01, "p < 0.01", sprintf("p = %.2f", p)))

norm_tool_label <- function(raw) {
  x <- as.character(raw)
  x <- gsub("_HMM[0-9]+", "", x); x <- gsub("_TumorFreeNormal", "", x)
  x <- gsub("_TumorFreeTME", "", x); x <- gsub("_TumorFree", "", x)
  x <- gsub("_Normal$", "", x); x <- gsub("_TME$", "", x)
  id <- tool_id(x); lab <- unname(ID2LABEL_T1[id])
  ifelse(is.na(lab), unname(tool_label(x)), lab)
}
metric_labeller <- as_labeller(function(x) { o <- unname(METRIC_LABELS[x]); ifelse(is.na(o), x, o) })

write_tab <- function(df, name) readr::write_csv(df, file.path(OUT_TAB, name))

# =====================================================================
# SYNTHETIC panels a / b   (input_data/synth/benchmark_metrics_normalref.csv)
# =====================================================================
raw <- read_csv(file.path(IN_SYN, "benchmark_metrics_normalref.csv"))
metric_cols   <- intersect(METRIC_ALL, colnames(raw))
rep_levels    <- sort(unique(raw$rep))
rep_label_map <- setNames(paste0("R", seq_along(rep_levels)), rep_levels)

dat <- raw %>%
  mutate(tool_lab = factor(norm_tool_label(.data$tool), levels = TOOL_ORDER_T1),
         scope    = factor(.data$scope, levels = SCOPE_LEVELS),
         rep_lab  = factor(unname(rep_label_map[as.character(.data$rep)]),
                           levels = unname(rep_label_map))) %>%
  select(tool_lab, scope, rep, rep_lab, all_of(metric_cols)) %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  filter(is.finite(.data$value)) %>%
  mutate(metric = factor(.data$metric, levels = metric_cols)) %>%
  filter(!is.na(.data$tool_lab))

METRIC_2A   <- c("spearman", "max_f1")
METRIC_MAIN <- c("max_f1")

sm_labels <- function(d) {
  d %>% distinct(.data$metric, .data$rep, .data$tool_lab, .data$value) %>%
    filter(is.finite(.data$value)) %>%
    group_by(.data$metric) %>%
    group_modify(~ tibble(p = skillings_mack(.x$value, as.character(.x$tool_lab),
                                             as.character(.x$rep))$p)) %>%
    ungroup() %>%
    transmute(.data$metric, label = fmt_p_disp(.data$p),
              tool_lab = factor(tail(TOOL_ORDER_T1, 1), levels = TOOL_ORDER_T1), y = Inf)
}

# ---- panel a ----
overall <- dat %>% filter(.data$scope == "overall", .data$metric %in% METRIC_2A) %>%
  mutate(metric = factor(.data$metric, levels = METRIC_2A))
sm_lab_a <- sm_labels(overall)
p_a <- ggplot(overall, aes(.data$tool_lab, .data$value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = .data$tool_lab), width = 0.15, height = 0, size = 0.55, alpha = 0.75) +
  scale_color_manual(values = TOOL_COLORS_T1, drop = FALSE, guide = "none") +
  geom_text(data = sm_lab_a, aes(.data$tool_lab, .data$y, label = .data$label),
            inherit.aes = FALSE, hjust = 1, vjust = 1.3, size = P_ANNO_SIZE,
            fontface = P_ANNO_FACE, family = FONT, color = "grey25", na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2, labeller = metric_labeller) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) +
  theme_pub(base_size = BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = BASE_SIZE - 1))
write_tab(overall %>% transmute(metric, tool_lab, rep, value), "fig2_panel_a_overall_accuracy.csv")
write_tab(sm_lab_a %>% transmute(metric, skillingsmack_label = label), "fig2_panel_a_skillingsmack.csv")

# ---- panel b ----
finish_box_v <- function(p, stat) {
  p +
    geom_text(data = stat, aes(.data$tool_lab, .data$y, label = .data$lab),
              inherit.aes = FALSE, hjust = 0.5, vjust = 1.2, size = P_ANNO_SIZE,
              fontface = P_ANNO_FACE, family = FONT, na.rm = TRUE) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(x = NULL, y = NULL) + theme_pub(base_size = BASE_SIZE) +
    theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(1, "pt"), legend.key.size = unit(0.20, "cm"),
          legend.text = element_text(size = BASE_SIZE - 2),
          strip.text = element_text(size = BASE_SIZE - 1, face = "bold"),
          axis.text.x = element_text(angle = 40, hjust = 1, size = BASE_SIZE - 2),
          axis.text.y = element_text(size = BASE_SIZE - 1))
}
scope_df <- dat %>% filter(.data$scope %in% c("broad_only", "narrow_only")) %>%
  filter(.data$metric %in% METRIC_MAIN) %>%
  mutate(metric = factor(as.character(.data$metric), levels = METRIC_MAIN)) %>% droplevels()
scope_cols <- c(broad_only = "#8DBEEA", narrow_only = "#EE9F9C")
scope_sig <- scope_df %>%
  select(metric, tool_lab, rep, scope, value) %>%
  pivot_wider(names_from = "scope", values_from = "value") %>%
  group_by(.data$metric, .data$tool_lab) %>%
  summarise(p = { ok <- is.finite(.data$broad_only) & is.finite(.data$narrow_only)
            if (sum(ok) < 2) NA_real_ else tryCatch(
              wilcox.test(.data$broad_only[ok], .data$narrow_only[ok], paired = TRUE,
                          exact = FALSE)$p.value, error = function(e) NA_real_) },
            .groups = "drop") %>%
  group_by(.data$metric) %>%
  mutate(p_adj = p.adjust(.data$p, method = "BH"),
         lab = dplyr::case_when(is.na(.data$p_adj) ~ "", .data$p_adj <= 0.001 ~ "***",
                                .data$p_adj <= 0.01 ~ "**", .data$p_adj <= 0.05 ~ "*",
                                TRUE ~ "ns")) %>%
  ungroup() %>% mutate(y = Inf)
p_b <- finish_box_v(
  ggplot(scope_df, aes(.data$tool_lab, .data$value, fill = .data$scope)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, outlier.shape = NA,
                 color = "black", linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               size = 0.4, alpha = 0.35, color = "black") +
    scale_fill_manual(values = scope_cols,
                      labels = c(broad_only = "Broad", narrow_only = "Narrow"), name = NULL),
  scope_sig) + labs(title = "Broad vs narrow CNVs", y = "Max F1") +
  theme(plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0.5, margin = margin(b = 1)),
        axis.title.y = element_text(size = BASE_SIZE - 1, margin = margin(r = 1)))
write_tab(scope_df %>% transmute(metric, tool_lab, rep, scope, value), "fig2_panel_b_broad_vs_narrow.csv")
write_tab(scope_sig %>% transmute(tool_lab, p_broad_vs_narrow = p, p_adj_BH = p_adj, star = lab),
          "fig2_panel_b_broad_vs_narrow_wilcox.csv")

# =====================================================================
# panel c : per-tool CNV accuracy vs tumor purity (Max F1)  [synth purity]
# =====================================================================
pur_raw2regid <- c(
  "InferCNV_HMM17_TumorFreeNormal_CNV"="infercnv_cnv","InferCNV_HMM17_TumorFreeNormal_Expr"="infercnv_expr",
  "InferCNVPy_TumorFreeNormal"="infercnvpy","Numbat_TumorFreeNormal_CNV"="numbat_cnv",
  "Numbat_TumorFreeNormal_Expr"="numbat_expr","CONICSmat_TumorFreeNormal"="conicsmat",
  "CopyKAT_TumorFreeNormal"="copykat","SCEVAN_TumorFreeNormal"="scevan",
  "XClone_TumorFreeNormal"="xclone","CalicoST"="calicost","CopyVAE"="copyvae","STARCH"="starch")
pur_raw2lab <- c(
  "InferCNV_HMM17_TumorFreeNormal_CNV"="InferCNV-CNV","InferCNV_HMM17_TumorFreeNormal_Expr"="InferCNV-Expr",
  "InferCNVPy_TumorFreeNormal"="InferCNVpy","Numbat_TumorFreeNormal_CNV"="Numbat-CNV",
  "Numbat_TumorFreeNormal_Expr"="Numbat-Expr","CONICSmat_TumorFreeNormal"="CONICSmat",
  "CopyKAT_TumorFreeNormal"="CopyKAT","SCEVAN_TumorFreeNormal"="SCEVAN",
  "XClone_TumorFreeNormal"="XClone","CalicoST"="CalicoST","CopyVAE"="CopyVAE","STARCH"="STARCH")
pur_col  <- setNames(t1$color, t1$id)
PUR_HI   <- c("InferCNV-CNV", "Numbat-CNV", "Numbat-Expr", "CalicoST")
PUR_GREY <- "#bdbdbd"
pur_dat <- read_csv(file.path(IN_SYN, "purity_bin_metrics_all_tools.csv")) %>%
  transmute(tool_raw = .data$tool, regid = unname(pur_raw2regid[.data$tool_raw]),
            lab = unname(pur_raw2lab[.data$tool_raw]), purity_bin = .data$purity_bin,
            purity_mid = as.numeric(sub("-.*$", "", .data$purity_bin)) + 0.05,
            spearman = .data$spearman, max_f1 = .data$max_f1) %>%
  pivot_longer(c("spearman", "max_f1"), names_to = "metric", values_to = "acc") %>%
  mutate(is_hi = .data$lab %in% PUR_HI,
         col = ifelse(.data$is_hi, unname(pur_col[.data$regid]), PUR_GREY))
pur_med <- pur_dat %>%
  group_by(.data$metric, .data$purity_bin, .data$purity_mid) %>%
  summarise(acc = median(.data$acc), .groups = "drop") %>%
  mutate(lab = "median", col = "black")
PUR_BINLABS <- pur_dat %>% distinct(.data$purity_mid, .data$purity_bin) %>% arrange(.data$purity_mid)

pur_panel <- function(mtr, ylo, yhi, ylab) {
  dd <- pur_dat %>% filter(.data$metric == mtr); mm <- pur_med %>% filter(.data$metric == mtr)
  others <- dd %>% filter(!.data$is_hi); hi <- dd %>% filter(.data$is_hi)
  ends <- bind_rows(hi %>% group_by(.data$lab) %>% slice_max(.data$purity_mid, n = 1) %>% ungroup(),
                    mm %>% slice_max(.data$purity_mid, n = 1))
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.05) +
    geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, color = "grey40") +
    geom_line(data = others, aes(.data$purity_mid, .data$acc, group = .data$lab, color = .data$col),
              linewidth = 0.4, alpha = 0.55) +
    geom_line(data = mm, aes(.data$purity_mid, .data$acc, color = .data$col),
              linetype = "dashed", linewidth = 0.9) +
    geom_line(data = hi, aes(.data$purity_mid, .data$acc, group = .data$lab, color = .data$col),
              linewidth = 0.9, alpha = 0.6) +
    geom_point(data = hi, aes(.data$purity_mid, .data$acc, color = .data$col),
               size = 1.0, alpha = 0.6, show.legend = FALSE) +
    geom_text_repel(data = ends, aes(.data$purity_mid, .data$acc, label = .data$lab, color = .data$col),
                    family = FONT, hjust = 0, direction = "y", nudge_x = 0.03, size = 2.0, segment.size = 0.2,
                    min.segment.length = 0, xlim = c(0.99, NA), show.legend = FALSE) +
    scale_color_identity(guide = "legend", name = NULL, breaks = c("black", PUR_GREY),
                         labels = c("median (all 12 tools)", "other tools (n=8)")) +
    scale_x_continuous(breaks = PUR_BINLABS$purity_mid, labels = PUR_BINLABS$purity_bin,
                       expand = expansion(mult = c(0.03, 0.34))) +
    coord_cartesian(ylim = c(ylo, yhi)) + labs(x = "Tumor proportion (bin)", y = ylab) +
    theme_pub(base_size = BASE_SIZE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
          legend.position = "bottom")
}
p_c <- pur_panel("max_f1", 0.33, 0.90, "Max F1") +
  theme(legend.position = "right", legend.text = element_text(size = BASE_SIZE - 2, family = FONT))
write_tab(pur_dat %>% filter(.data$metric == "max_f1") %>%
            transmute(tool = lab, regid, purity_bin, purity_mid, max_f1 = acc, highlighted = is_hi),
          "fig2_panel_c_purity_curve_maxf1.csv")

# =====================================================================
# REAL panels d / e / f    (input_data/real/pseudobulk_vs_wes_metrics_all_datasets.csv via final_record_select_task1())
# BEST.csv is the frozen best-run-per-(sample,tool) table = the `best0` the
# original computed from all_datasets.csv, so it is consumed directly here.
# =====================================================================
METRICS4    <- c("Spearman", "max_F1", "AUC_gain_trunc", "AUC_loss_trunc")
METRICS_2EF <- c("Spearman", "max_F1")
METRIC_2F   <- "max_F1"
RMETRIC_LABELS <- c(Spearman = "Correlation", max_F1 = "Max F1",
                    AUC_gain_trunc = "AUC gain", AUC_loss_trunc = "AUC loss")
rmetric_labeller <- as_labeller(function(x) { o <- unname(RMETRIC_LABELS[x]); ifelse(is.na(o), x, o) })
CNV_NEGATIVE <- c("skin_P4_V1", "skin_P4_V2")

ds <- figure_settings("dataset", "global")
DSID2LABEL    <- setNames(ds$label, ds$id)
DATASET_ORDER <- unname(DSID2LABEL[c("Ding_nature_wes", "metastatic_prostate_ffpe", "skin_tumor_ff")])
DING_LABEL    <- unname(DSID2LABEL["Ding_nature_wes"])
.alias_lc <- tibble(key = tolower(TOOL_ALIAS_TABLE$alias), id = TOOL_ALIAS_TABLE$id) %>%
  distinct(key, .keep_all = TRUE)
norm_tool_id <- function(x) { id <- .alias_lc$id[match(tolower(as.character(x)), .alias_lc$key)]
  ifelse(is.na(id), as.character(x), id) }

best0 <- final_record_select_task1(read_csv(file.path(IN_REAL, "pseudobulk_vs_wes_metrics_all_datasets.csv"))) %>%
  mutate(sample = trimws(as.character(.data$sample)),
         tool_id = norm_tool_id(.data$tool),
         tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab),
         dataset_lab = unname(DSID2LABEL[as.character(.data$dataset)]),
         sample_type = ifelse(tolower(trimws(.data$sample_type)) == "ffpe", "ffpe", "ff")) %>%
  filter(!is.na(.data$dataset_lab), !.data$sample %in% CNV_NEGATIVE)
# BEST.csv should already be one row per (dataset, sample, tool_id):
stopifnot(!any(duplicated(best0[c("dataset", "sample", "tool_id")])))

present_meta  <- TOOL_TABLE %>% filter(.data$id %in% unique(best0$tool_id)) %>% arrange(.data$order)
TOOL_ORDER_R  <- present_meta$label
TOOL_COLORS_R <- setNames(present_meta$color, present_meta$label)
best0 <- best0 %>% mutate(tool_lab = factor(.data$tool_lab, levels = TOOL_ORDER_R),
                          dataset_lab = factor(.data$dataset_lab, levels = DATASET_ORDER))

# ---- panel d : accuracy across datasets (Max F1) ----
f_long <- best0 %>%
  pivot_longer(all_of(METRIC_2F), names_to = "metric", values_to = "value") %>%
  filter(is.finite(.data$value)) %>% mutate(metric = factor(.data$metric, levels = METRIC_2F))
f_levels <- intersect(levels(f_long$dataset_lab), unique(f_long$dataset_lab))
f_comparisons <- if (length(f_levels) >= 2) utils::combn(f_levels, 2, simplify = FALSE) else list()
f_ymax <- max(f_long$value, na.rm = TRUE); f_pad <- 0.05 * (f_ymax - min(f_long$value, na.rm = TRUE))
p_d <- ggplot(f_long, aes(.data$dataset_lab, .data$value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = .data$tool_lab), width = 0.18, height = 0, size = 0.9, alpha = 0.85) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y", labeller = rmetric_labeller) +
  scale_color_manual(values = TOOL_COLORS_R, drop = FALSE, guide = "none") +
  labs(title = "Accuracy across datasets", x = NULL, y = NULL) +
  theme_pub(base_size = BASE_SIZE) +
  theme(plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0.5, margin = margin(b = 1)),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1, size = BASE_SIZE - 1))
f_brk <- NULL
if (length(f_comparisons) > 0) {
  f_by_tool <- f_long %>% group_by(.data$dataset_lab, .data$tool_lab) %>%
    summarise(v = median(.data$value, na.rm = TRUE), .groups = "drop")
  f_rng <- f_ymax - min(f_long$value, na.rm = TRUE)
  f_pvals <- vapply(f_comparisons, function(pr) {
    w <- f_by_tool %>% filter(.data$dataset_lab %in% pr) %>%
      tidyr::pivot_wider(names_from = "dataset_lab", values_from = "v")
    a <- w[[pr[[1]]]]; b <- w[[pr[[2]]]]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 2L) return(NA_real_)
    tryCatch(stats::wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  }, numeric(1))
  f_ypos <- f_ymax + f_rng * (0.07 + 0.15 * (seq_along(f_comparisons) - 1)); f_tick <- f_rng * 0.02
  f_brk <- data.frame(
    x = vapply(f_comparisons, function(pr) match(pr[[1]], f_levels), numeric(1)),
    xend = vapply(f_comparisons, function(pr) match(pr[[2]], f_levels), numeric(1)),
    y = f_ypos, lab = vapply(f_pvals, fmt_p_disp, character(1)),
    metric = factor(METRIC_2F, levels = METRIC_2F))
  p_d <- p_d +
    geom_segment(data = f_brk, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
                 inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
    geom_segment(data = f_brk, aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$y - f_tick),
                 inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
    geom_segment(data = f_brk, aes(x = .data$xend, xend = .data$xend, y = .data$y, yend = .data$y - f_tick),
                 inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
    geom_text(data = f_brk, aes(x = (.data$x + .data$xend) / 2, y = .data$y, label = .data$lab),
              inherit.aes = FALSE, vjust = -0.35, size = P_ANNO_SIZE, family = FONT) +
    expand_limits(y = max(f_ypos) + f_rng * 0.14)
}
write_tab(f_long %>% transmute(dataset_lab, sample, tool_lab, metric, value),
          "fig2_panel_d_accuracy_across_datasets.csv")
if (!is.null(f_brk)) write_tab(
  data.frame(dataset1 = vapply(f_comparisons, `[[`, character(1), 1),
             dataset2 = vapply(f_comparisons, `[[`, character(1), 2),
             p_paired_wilcox_by_tool = f_pvals),
  "fig2_panel_d_between_dataset_wilcox.csv")

# ---- panels e / f : Ding FFPE / Ding FF, tools on y, 2 metrics ----
# inline block_omnibus (Friedman if the tool x sample grid is complete, else
# Skillings-Mack), matching the original fig2 shared stats lib.
.ginv_svd <- function(M, tol = 1e-9) {
  s <- svd(M); pos <- s$d > tol * max(s$d)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}
sm_matrix <- function(mat) {
  obs <- !is.na(mat); keep <- rowSums(obs) >= 2L
  mat <- mat[keep, , drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  g <- ncol(mat); b <- nrow(mat)
  A <- setNames(numeric(g), colnames(mat)); Sigma <- matrix(0, g, g, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(b)) {
    present <- which(obs[j, ]); k <- length(present)
    r <- rank(mat[j, present], ties.method = "average")
    a <- sqrt(12 / (k + 1)) * (r - (k + 1) / 2)
    A[present] <- A[present] + a
    Sigma[present, present] <- Sigma[present, present] - 1
    diag(Sigma)[present] <- diag(Sigma)[present] + k
  }
  stat <- as.numeric(t(A) %*% .ginv_svd(Sigma) %*% A); df <- g - 1L
  list(method = "Skillings-Mack", p.value = pchisq(stat, df, lower.tail = FALSE))
}
block_omnibus <- function(d) {
  d <- as.data.frame(d); d <- d[!is.na(d$block) & !is.na(d$group), c("value", "group", "block")]
  mat <- tapply(d$value, list(d$block, d$group), function(x) mean(x, na.rm = TRUE)); mat[is.nan(mat)] <- NA
  if (all(!is.na(mat))) {
    dd <- d; dd$group <- factor(dd$group); dd$block <- factor(dd$block)
    ft <- friedman.test(value ~ group | block, data = dd)
    return(list(method = "Friedman", p.value = ft$p.value))
  }
  sm_matrix(mat)
}
OMNI_METRIC <- c(Spearman = "spearman", max_F1 = "maxf1")

make_tool_box <- function(df, metrics, panel_title, order_metric = "max_F1") {
  keep <- intersect(METRICS4, metrics)
  d <- df %>% pivot_longer(all_of(keep), names_to = "metric", values_to = "value") %>%
    filter(is.finite(.data$value)) %>%
    mutate(metric = factor(.data$metric, levels = keep),
           tool_lab = droplevels(factor(.data$tool_lab, levels = TOOL_ORDER_R)))
  med <- d %>% filter(.data$metric == order_metric) %>%
    group_by(.data$tool_lab) %>% summarise(m = median(.data$value, na.rm = TRUE), .groups = "drop") %>%
    arrange(dplyr::desc(.data$m))
  ord <- c(as.character(med$tool_lab), setdiff(levels(d$tool_lab), as.character(med$tool_lab)))
  d <- d %>% mutate(tool_lab = factor(as.character(.data$tool_lab), levels = ord))
  # omnibus p per metric
  omni <- do.call(rbind, lapply(keep, function(mc) {
    dd <- data.frame(value = df[[mc]], group = as.character(df$tool_lab), block = as.character(df$sample))
    dd <- dd[is.finite(dd$value), , drop = FALSE]
    if (length(unique(dd$group)) < 2 || nrow(dd) < 2) return(NULL)
    data.frame(metric = mc, p = block_omnibus(dd)$p.value)
  }))
  anno <- NULL
  if (!is.null(omni)) anno <- omni %>% transmute(metric = factor(metric, levels = keep),
                                                 label = fmt_p_disp(p))
  p <- ggplot(d, aes(.data$tool_lab, .data$value, fill = .data$tool_lab)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, color = "grey25", linewidth = 0.35) +
    facet_wrap(~ metric, ncol = length(keep), scales = "free_y", labeller = rmetric_labeller) +
    scale_fill_manual(values = TOOL_COLORS_R, drop = FALSE, guide = "none") +
    scale_y_continuous(breaks = scales::extended_breaks(n = 3)) +
    coord_flip(clip = "off") +
    scale_x_discrete(limits = rev, expand = expansion(add = c(0.6, 1.5))) +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_pub(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0.5, margin = margin(b = 1)),
          strip.text = element_text(size = BASE_SIZE - 2, face = "bold"),
          axis.text.x = element_text(size = BASE_SIZE - 2), axis.text.y = element_text(size = BASE_SIZE - 2),
          panel.spacing = unit(2, "pt"),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3))
  if (!is.null(anno))
    p <- p + geom_text(data = anno, aes(x = Inf, y = Inf, label = .data$label),
                       hjust = 1.05, vjust = 1.1, size = P_ANNO_SIZE, fontface = P_ANNO_FACE,
                       family = FONT, inherit.aes = FALSE)
  list(plot = p, data = d, omni = omni)
}
ding_ffpe <- best0 %>% filter(.data$dataset_lab == DING_LABEL, .data$sample_type == "ffpe")
ding_ff   <- best0 %>% filter(.data$dataset_lab == DING_LABEL, .data$sample_type == "ff")
res_e <- make_tool_box(ding_ffpe, METRICS_2EF, paste0(DING_LABEL, " — FFPE"))
res_f <- make_tool_box(ding_ff,   METRICS_2EF, paste0(DING_LABEL, " — FF (expr only)"))
p_e <- res_e$plot; p_f <- res_f$plot
write_tab(res_e$data %>% transmute(sample, tool_lab, metric, value), "fig2_panel_e_ding_ffpe.csv")
write_tab(res_f$data %>% transmute(sample, tool_lab, metric, value), "fig2_panel_f_ding_ff.csv")
if (!is.null(res_e$omni)) write_tab(res_e$omni, "fig2_panel_e_ding_ffpe_omnibus.csv")
if (!is.null(res_f$omni)) write_tab(res_f$omni, "fig2_panel_f_ding_ff_omnibus.csv")

# Panel-e (Ding FFPE) ranking best->worst by median Max F1 (reused by heatmap rows).
G_RANK_IDS <- ding_ffpe %>% group_by(.data$tool_id) %>%
  summarise(m = median(.data$max_F1, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(.data$m)) %>% arrange(dplyr::desc(.data$m)) %>% pull(.data$tool_id)

# ---- shared real-data tool legend (keys d; identifies e/f/g tools) ----
leg_df <- tibble(tool_lab = factor(TOOL_ORDER_R, levels = TOOL_ORDER_R), x = 1, y = seq_along(TOOL_ORDER_R))
p_leg_full <- ggplot(leg_df, aes(.data$x, .data$y, color = .data$tool_lab)) +
  geom_point(size = 1.4) +
  scale_color_manual(values = TOOL_COLORS_R, drop = FALSE, name = "Tool") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.8))) +
  theme_pub(base_size = BASE_SIZE) +
  theme(legend.position = "right",
        legend.title = element_text(face = "plain", size = BASE_SIZE - 1, family = FONT),
        legend.text = element_text(size = BASE_SIZE - 2, family = FONT),
        legend.key.size = unit(0.28, "cm"))
tool_legend <- {
  g <- if (requireNamespace("cowplot", quietly = TRUE)) cowplot::get_legend(p_leg_full) else NULL
  if (is.null(g)) patchwork::guide_area() else patchwork::wrap_elements(full = g)
}

# =====================================================================
# panel g : gene-level CNV heatmap, HT397B1-H2U1
# tool rows from input_data/real/pseudobulk_gene/HT397B1-H2U1/<tool_id>.csv
# =====================================================================
STATE_TOOL_IDS <- c("infercnv_cnv", "numbat_cnv", "starch"); CN_LIMS <- c(1, 3)
CHR_ORDER <- paste0("chr", 1:22)
SAMPLE <- "HT397B1-H2U1"; TRUTH_WES_ID <- "merged_HT397B1-S1H4A4D1_1-12"
PSEUDO_DIR <- file.path(IN_REAL, "pseudobulk_gene", SAMPLE)
log2ratio_to_cn <- function(x) 2 * (2 ^ x); state_to_cn <- function(x) 2 + x
normalize_chr <- function(x) { x <- as.character(x); x <- sub("^CHR", "chr", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE); paste0("chr", x) }

gene_tbl <- read.delim(file.path(IN_REAL, "hgTables_hg38_gencode.txt"), sep = "\t",
                       header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_tbl <- gene_tbl[, c("name2", "chrom", "cdsStart", "cdsEnd")]
gene_tbl$mid <- (gene_tbl$cdsStart + gene_tbl$cdsEnd) %/% 2L
gene_tbl <- gene_tbl[gene_tbl$chrom %in% CHR_ORDER, ] %>%
  group_by(.data$name2) %>% slice_max(.data$cdsEnd - .data$cdsStart, n = 1, with_ties = FALSE) %>%
  ungroup() %>% mutate(chrom = factor(.data$chrom, levels = CHR_ORDER)) %>%
  arrange(.data$chrom, .data$mid) %>% mutate(gene_order = row_number())
chr_boundaries <- gene_tbl %>% group_by(.data$chrom) %>%
  summarise(first = min(.data$gene_order), last = max(.data$gene_order), .groups = "drop") %>%
  mutate(mid = (.data$first + .data$last) / 2)
chr_vlines <- chr_boundaries %>% mutate(end = .data$last + 0.5) %>%
  filter(.data$chrom != tail(CHR_ORDER, 1)) %>% pull(.data$end)

load_pseudo <- function(path, tid) {
  d <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- "gene"; val_col <- setdiff(names(d), "gene")[1]
  raw_v <- suppressWarnings(as.numeric(d[[val_col]]))
  cn <- if (tid %in% STATE_TOOL_IDS) state_to_cn(raw_v) else log2ratio_to_cn(raw_v)
  data.frame(gene = as.character(d$gene), copy_number = cn, tool_id = tid, stringsAsFactors = FALSE)
}
tool_files <- list.files(PSEUDO_DIR, pattern = "\\.csv$", full.names = TRUE)
tool_long <- purrr::map_dfr(tool_files, function(p) load_pseudo(p, sub("\\.csv$", "", basename(p))))

seg <- read.csv(file.path(IN_REAL, "HT397B1-H2U1_wes_cnv.csv"), header = TRUE,
                check.names = FALSE, stringsAsFactors = FALSE)
seg <- seg[as.character(seg$sample) == TRUTH_WES_ID, , drop = FALSE]
seg$CONTIG <- normalize_chr(seg$CONTIG); seg <- seg[seg$CONTIG %in% CHR_ORDER, ]
wes <- gene_tbl %>% select(.data$name2, .data$chrom, .data$mid, .data$gene_order); wes$log2ratio <- NA_real_
for (chr in CHR_ORDER) {
  g_idx <- which(wes$chrom == chr); if (!length(g_idx)) next
  s <- seg[seg$CONTIG == chr, ]; if (!nrow(s)) next
  s <- s[order(s$START), ]
  starts <- as.integer(round(as.numeric(s$START))); ends <- as.integer(round(as.numeric(s$END)))
  vals <- suppressWarnings(as.numeric(s$MEAN_LOG2_COPY_RATIO)); g_mid <- wes$mid[g_idx]
  idx <- findInterval(g_mid, starts); hit <- idx >= 1 & idx <= length(ends)
  in_seg <- hit & (g_mid <= ends[pmax(1, idx)]); wes$log2ratio[g_idx[in_seg]] <- vals[idx[in_seg]]
}
wes_long <- data.frame(gene = wes$name2, copy_number = log2ratio_to_cn(wes$log2ratio),
                       tool_id = "WES", stringsAsFactors = FALSE)
all_long <- bind_rows(wes_long, tool_long) %>%
  inner_join(gene_tbl %>% select(.data$name2, .data$gene_order, .data$chrom), by = c("gene" = "name2"))
tool_labels <- setNames(TOOL_TABLE$label, TOOL_TABLE$id)

tool_ids_present <- setdiff(unique(all_long$tool_id), "WES")
ranked <- intersect(G_RANK_IDS, tool_ids_present)
rest   <- TOOL_TABLE %>% filter(.data$id %in% setdiff(tool_ids_present, ranked)) %>%
  arrange(.data$order) %>% pull(.data$id)
tool_ids_used <- c(ranked, rest)
row_levels <- c("WES", tool_ids_used)
row_labels <- setNames(c("WES (truth)", unname(tool_labels[tool_ids_used])), row_levels)
all_long <- all_long %>% mutate(row_label = unname(row_labels[.data$tool_id]),
                                row_id = factor(.data$row_label, levels = rev(row_labels)))
p_g <- ggplot(all_long, aes(.data$gene_order, .data$row_id, fill = .data$copy_number)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = chr_vlines, color = "grey55", linewidth = 0.18) +
  scale_fill_gradient2(low = "#2B5CB8", mid = "white", high = "#C9342C", midpoint = 2,
                       limits = CN_LIMS, breaks = c(1, 2, 3), labels = c("Loss", "Neutral", "Gain"),
                       oob = scales::squish, na.value = "grey92", name = "CNV\nstate") +
  scale_x_continuous(breaks = chr_boundaries$mid,
                     labels = sub("^chr", "", as.character(chr_boundaries$chrom)), expand = c(0, 0)) +
  labs(title = paste0("HT397B1-H2U1 — ", DING_LABEL, ", same-block WES"), x = "Chromosome", y = NULL) +
  theme_pub(base_size = BASE_SIZE) +
  theme(plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0, margin = margin(b = 1)),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = BASE_SIZE - 2), axis.text.y = element_text(size = BASE_SIZE - 1),
        legend.position = "right", legend.key.size = unit(0.30, "cm"),
        legend.title = element_text(size = BASE_SIZE - 1), legend.text = element_text(size = BASE_SIZE - 1))
write_tab(all_long %>% transmute(row_label, tool_id, gene, chrom, gene_order, copy_number),
          "fig2_panel_g_heatmap_HT397B1-H2U1.csv")

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
TAG <- theme(plot.tag = element_text(face = "bold", size = 8))
fr <- function(x) patchwork::free(x, type = "space", side = "l")
row_b  <- p_b + labs(tag = "b")
row_cd <- (p_c + labs(tag = "c")) | (p_d + labs(tag = "d")) | tool_legend
row_cd <- row_cd + plot_layout(widths = c(1.45, 0.78, 0.42))
row_ef <- (p_e + labs(tag = "e")) | (p_f + labs(tag = "f")); row_ef <- row_ef + plot_layout(widths = c(1, 1))
row_g  <- p_g + labs(tag = "g")
final <- (p_a + labs(tag = "a")) / fr(row_b) / fr(row_cd) / row_ef / row_g +
  plot_layout(heights = c(1.05, 1.6, 1.75, 1.6, 1.4)) &
  TAG & theme(text = element_text(family = FONT),
              axis.text.x = element_text(margin = margin(t = 0.5)),
              axis.text.y = element_text(margin = margin(r = 0.5)))

W <- 7.09; H <- 8.85
ggsave(file.path(OUT_FIG, "fig2_compose.pdf"), final, width = W, height = H, limitsize = FALSE, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig2_compose.png"), final, width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)

finish <- theme(text = element_text(family = FONT),
                axis.text.x = element_text(margin = margin(t = 0.5)),
                axis.text.y = element_text(margin = margin(r = 0.5)))
panels <- list(a = list(p_a, 7.09, 2.3), b = list(p_b, 7.09, 2.6), c = list(p_c, 3.9, 2.9),
               d = list(p_d, 3.2, 2.9), e = list(p_e, 3.6, 3.2), f = list(p_f, 3.6, 3.2),
               g = list(p_g, 7.09, 2.6))
for (nm in names(panels)) {
  pl <- panels[[nm]]; p <- pl[[1]] + finish; base <- file.path(OUT_FIG, paste0("fig2_panel_", nm))
  ggsave(paste0(base, ".pdf"), p, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] fig2_compose + panels a-g (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
