#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig S_1A (Task 1, SYNTHETIC diagnostics) — rebuilt from scratch.
#
# Reproduces paper/figures/supp_task1_merge/figS_task1A_syn_overall_compose.R
# using ONLY input_data/synth. Panels:
#   a Overall accuracy by tool type (2x2 metric grid), group-wise omnibus p
#   b Pairwise paired-Wilcoxon heatmaps (Spearman, Max F1, pAUC gain/loss)
#   c Broad vs narrow CNV recovery (Spearman, pAUC gain/loss), 3x1 strip
#   d Per-tool accuracy vs tumor purity (Spearman)
#
# Sources only utils/figure_setup.R (settings + theme + skillings_mack).
# Outputs: composed fig + panels (pdf/png) in supp_figures/, stats + plotted
# data in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(readr); library(purrr); library(ggrepel); library(tibble)
})

.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]))); normalizePath("figS_task1A_syn.R")
})()
SCRIPTS_DIR <- dirname(.this_file); TASK_DIR <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS   <- normalizePath(file.path(TASK_DIR, "..", "utils", "figure_setup.R"))
IN_SYN  <- file.path(TASK_DIR, "input_data", "synth")
OUT_FIG <- file.path(TASK_DIR, "supp_figures"); OUT_TAB <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(UTILS)
write_tab <- function(df, name) readr::write_csv(df, file.path(OUT_TAB, name))

FONT <- pick_font()   # real "Arial" if installed, else bundled "Liberation Sans" (Arial metric clone), else "sans"
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

METRIC_ORDER <- c("spearman", "max_f1", "pauc_gain", "pauc_loss")
SUPP_METRICS <- c("spearman", "pauc_gain", "pauc_loss")
SCOPE_LEVELS <- c("overall", "broad_only", "narrow_only")
BASE <- 7; P_ANNO_SIZE <- 5.5 / .pt; P_ANNO_FACE <- "plain"

t1 <- figure_settings("tool_order", "synth_task1")
TOOL_ORDER_T1 <- t1$label; TOOL_COLORS_T1 <- setNames(t1$color, t1$label)
ID2LABEL_T1 <- setNames(t1$label, t1$id)
METRIC_LABELS <- setting_vector("metric", "global", "label"); METRIC_LABELS["spearman"] <- "Spearman"
metric_labeller <- as_labeller(function(x) { o <- unname(METRIC_LABELS[as.character(x)]); ifelse(is.na(o), as.character(x), o) })

norm_tool_label <- function(raw) {
  x <- as.character(raw)
  x <- gsub("_HMM[0-9]+", "", x); x <- gsub("_TumorFreeNormal", "", x); x <- gsub("_TumorFreeTME", "", x)
  x <- gsub("_TumorFree", "", x); x <- gsub("_Normal$", "", x); x <- gsub("_TME$", "", x)
  id <- tool_id(x); lab <- unname(ID2LABEL_T1[id]); ifelse(is.na(lab), unname(tool_label(x)), lab)
}
star_labels <- function(p) dplyr::case_when(is.na(p) ~ "", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "ns")
sig_label   <- function(p) dplyr::case_when(is.na(p) ~ "", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "")
fmt_p_disp  <- function(p) dplyr::case_when(is.na(p) ~ "p = NA", p < 0.01 ~ "p < 0.01", TRUE ~ sprintf("p = %.2f", p))

raw <- read_csv(file.path(IN_SYN, "benchmark_metrics_normalref.csv"))
stopifnot(all(c("tool", "scope", "rep") %in% colnames(raw)))
dat <- raw %>%
  mutate(tool_lab = factor(norm_tool_label(.data$tool), levels = TOOL_ORDER_T1),
         scope = factor(.data$scope, levels = SCOPE_LEVELS)) %>%
  filter(!is.na(.data$tool_lab)) %>%
  select(tool_lab, scope, rep, all_of(METRIC_ORDER)) %>%
  pivot_longer(all_of(METRIC_ORDER), names_to = "metric", values_to = "value") %>%
  filter(is.finite(.data$value)) %>%
  mutate(metric = factor(.data$metric, levels = METRIC_ORDER))

# ---- (a) overall accuracy by tool type ------------------------------
make_block_matrix <- function(d, tools) {
  blocks <- sort(unique(as.character(d$rep)))
  M <- matrix(NA_real_, length(blocks), length(tools), dimnames = list(blocks, tools))
  d <- d %>% distinct(.data$rep, .data$tool_lab, .data$value)
  for (i in seq_len(nrow(d))) M[as.character(d$rep[i]), as.character(d$tool_lab[i])] <- d$value[i]
  M
}
group_omnibus <- function(d, tool_order) {
  d %>% distinct(.data$metric, .data$tool_type, .data$rep, .data$tool_lab, .data$value) %>%
    group_by(.data$metric, .data$tool_type) %>%
    group_modify(~ {
      tools <- tool_order[tool_order %in% as.character(unique(.x$tool_lab))]
      M <- make_block_matrix(.x, tools); n_missing <- sum(!is.finite(M))
      method <- if (n_missing == 0) "Friedman" else "Skillings-Mack"
      if (nrow(M) < 2 || ncol(M) < 2)
        return(tibble(method = method, statistic = NA_real_, df = NA_integer_, p = NA_real_))
      if (method == "Friedman") {
        ft <- tryCatch(stats::friedman.test(M), error = function(e) NULL)
        if (is.null(ft)) return(tibble(method = method, statistic = NA_real_, df = NA_integer_, p = NA_real_))
        return(tibble(method = method, statistic = unname(ft$statistic), df = unname(ft$parameter), p = ft$p.value))
      }
      sm <- skillings_mack(.x$value, as.character(.x$tool_lab), as.character(.x$rep))
      tibble(method = method, statistic = sm$statistic, df = sm$df, p = sm$p)
    }) %>% ungroup()
}
overall <- dat %>% filter(.data$scope == "overall") %>% droplevels()
ALLELE_TOOLS <- c("Numbat (CNV)", "XClone", "CalicoST")
observed_tools <- TOOL_ORDER_T1[TOOL_ORDER_T1 %in% as.character(unique(overall$tool_lab))]
expr_tools <- observed_tools[!observed_tools %in% ALLELE_TOOLS]
allele_tools <- observed_tools[observed_tools %in% ALLELE_TOOLS]
TOOL_ORDER_BY_TYPE <- c(expr_tools, allele_tools); divider_x <- length(expr_tools) + 0.5
overall_type <- overall %>%
  mutate(tool_lab = factor(as.character(.data$tool_lab), levels = TOOL_ORDER_BY_TYPE),
         tool_x = match(as.character(.data$tool_lab), TOOL_ORDER_BY_TYPE),
         tool_type = factor(if_else(as.character(.data$tool_lab) %in% ALLELE_TOOLS, "Allele-based", "Expression-based"),
                            levels = c("Expression-based", "Allele-based")))
type_omnibus <- group_omnibus(overall_type, TOOL_ORDER_BY_TYPE)
write_tab(type_omnibus %>% arrange(.data$metric, .data$tool_type), "figS_task1A_syn_type_omnibus.csv")
type_centers <- tibble(tool_type = factor(c("Expression-based", "Allele-based"), levels = c("Expression-based", "Allele-based")),
                       x = c((length(expr_tools) + 1) / 2, length(expr_tools) + (length(allele_tools) + 1) / 2))
p_type_lab <- type_omnibus %>% left_join(type_centers, by = "tool_type") %>% mutate(label = fmt_p_disp(.data$p), y = Inf)
pA_type <- ggplot(overall_type, aes(.data$tool_x, .data$value, group = .data$tool_lab)) +
  geom_boxplot(width = 0.68, outlier.shape = NA, fill = "white", color = "grey25", linewidth = 0.35) +
  geom_jitter(aes(color = .data$tool_lab), width = 0.13, height = 0, size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = divider_x, color = "grey25", linewidth = 0.35) +
  geom_text(data = p_type_lab, aes(.data$x, .data$y, label = .data$label), inherit.aes = FALSE,
            vjust = 1.15, size = P_ANNO_SIZE, lineheight = 0.85, fontface = P_ANNO_FACE, color = "grey25", na.rm = TRUE) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y", labeller = labeller(metric = metric_labeller)) +
  scale_x_continuous(breaks = seq_along(TOOL_ORDER_BY_TYPE), labels = TOOL_ORDER_BY_TYPE, expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  scale_color_manual(values = TOOL_COLORS_T1, drop = FALSE, guide = "none") +
  labs(title = "Synthetic accuracy by tool type", subtitle = "Expression-based (left) vs allele-based (right)", x = NULL, y = "Metric value") +
  theme_pub(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE - 2), axis.text.y = element_text(size = BASE - 1),
        strip.text = element_text(size = BASE - 1, face = "bold"), plot.title = element_text(size = BASE, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = BASE - 1, hjust = 0), panel.spacing = unit(2.5, "pt"))
write_tab(overall_type %>% transmute(tool_lab, tool_type, rep, metric, value), "figS_task1A_panel_a_accuracy_by_type.csv")

# ---- (b in layout / drawn as c-code) broad vs narrow ----------------
scope_df <- dat %>% filter(.data$scope %in% c("broad_only", "narrow_only"), .data$metric %in% SUPP_METRICS) %>%
  mutate(tool_lab = factor(as.character(.data$tool_lab), levels = TOOL_ORDER_T1),
         metric = factor(as.character(.data$metric), levels = SUPP_METRICS)) %>% droplevels()
scope_sig <- scope_df %>% select(metric, tool_lab, rep, scope, value) %>%
  pivot_wider(names_from = "scope", values_from = "value") %>%
  group_by(.data$metric, .data$tool_lab) %>%
  summarise(p = { ok <- is.finite(.data$broad_only) & is.finite(.data$narrow_only)
    if (sum(ok) < 2) NA_real_ else tryCatch(wilcox.test(.data$broad_only[ok], .data$narrow_only[ok], paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_) },
    .groups = "drop") %>%
  group_by(.data$metric) %>% mutate(p_adj = p.adjust(.data$p, method = "BH"), lab = star_labels(.data$p_adj)) %>% ungroup() %>% mutate(y = Inf)
scope_cols <- c(broad_only = "#8DBEEA", narrow_only = "#EE9F9C")
pB_scope <- ggplot(scope_df, aes(.data$tool_lab, .data$value, fill = .data$scope)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.64, outlier.shape = NA, color = "black", linewidth = 0.35) +
  geom_point(position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.8), size = 0.38, alpha = 0.30, color = "black") +
  geom_text(data = scope_sig, aes(.data$tool_lab, .data$y, label = .data$lab), inherit.aes = FALSE, vjust = 1.15, size = P_ANNO_SIZE, fontface = P_ANNO_FACE, na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, labeller = metric_labeller) +
  scale_fill_manual(values = scope_cols, labels = c(broad_only = "Broad", narrow_only = "Narrow"), name = NULL) +
  coord_cartesian(clip = "off") + labs(title = "Broad vs narrow\nCNV recovery", x = NULL, y = "Metric value") +
  theme_pub(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 2),
        plot.title = element_text(size = BASE, face = "bold", hjust = 0.5), legend.position = "top",
        legend.key.size = unit(0.25, "cm"), strip.text = element_text(size = BASE - 1, face = "bold"))
write_tab(scope_df %>% transmute(metric, tool_lab, rep, scope, value), "figS_task1A_panel_c_broad_vs_narrow.csv")
write_tab(scope_sig %>% transmute(metric, tool_lab, p, p_adj_BH = p_adj, star = lab), "figS_task1A_panel_c_broad_vs_narrow_wilcox.csv")

# ---- (b) pairwise paired-Wilcoxon heatmaps (2x2) --------------------
pairwise_for_metric <- function(d) {
  tools <- TOOL_ORDER_T1[TOOL_ORDER_T1 %in% as.character(unique(d$tool_lab))]
  pairs <- utils::combn(tools, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    wide <- d %>% filter(as.character(.data$tool_lab) %in% pair) %>%
      mutate(tool_lab = as.character(.data$tool_lab)) %>% select(rep, tool_lab, value) %>%
      pivot_wider(names_from = tool_lab, values_from = value)
    x <- wide[[pair[1]]]; y <- wide[[pair[2]]]; ok <- is.finite(x) & is.finite(y); n_pairs <- sum(ok)
    p_value <- if (n_pairs < 2) NA_real_ else tryCatch(wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
    tibble(group1 = pair[1], group2 = pair[2], n_pairs = n_pairs, p_value = p_value)
  })
}
pairwise <- overall %>% distinct(.data$metric, .data$rep, .data$tool_lab, .data$value) %>%
  group_by(.data$metric) %>% group_modify(~ pairwise_for_metric(.x)) %>% ungroup() %>%
  group_by(.data$metric) %>% mutate(p_adj_BH = p.adjust(.data$p_value, method = "BH"), signif = sig_label(.data$p_adj_BH)) %>% ungroup()
write_tab(pairwise %>% arrange(.data$metric, .data$group1, .data$group2), "figS_task1A_syn_pairwise_wilcox.csv")
make_sym <- function(d, tools) {
  bind_rows(transmute(d, row_tool = .data$group1, col_tool = .data$group2, signif = .data$signif, p_adj_BH = .data$p_adj_BH),
            transmute(d, row_tool = .data$group2, col_tool = .data$group1, signif = .data$signif, p_adj_BH = .data$p_adj_BH),
            tibble(row_tool = tools, col_tool = tools, signif = "", p_adj_BH = NA_real_)) %>%
    mutate(col_tool = factor(.data$col_tool, levels = tools), row_tool = factor(.data$row_tool, levels = rev(tools)),
           neg_log10_p = if_else(is.finite(.data$p_adj_BH), -log10(pmax(.data$p_adj_BH, .Machine$double.xmin)), NA_real_))
}
make_pair_plot <- function(metric_id, title, fill_limits) {
  sub <- pairwise %>% filter(.data$metric == metric_id)
  tools <- TOOL_ORDER_T1[TOOL_ORDER_T1 %in% unique(c(sub$group1, sub$group2))]
  heat <- make_sym(sub, tools)
  ggplot(heat, aes(.data$col_tool, .data$row_tool, fill = .data$neg_log10_p)) +
    geom_tile(color = "grey88", linewidth = 0.22) +
    geom_text(aes(label = .data$signif), size = 1.5, vjust = 0.75, color = "grey15", na.rm = TRUE) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = expression(-log[10](p[adj])), na.value = "grey95", limits = fill_limits, oob = scales::squish) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) + coord_equal() + labs(title = title, x = NULL, y = NULL) +
    theme_pub(base_size = BASE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 2), axis.text.y = element_text(size = BASE - 2),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.22, "cm"), legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
          plot.title = element_text(size = BASE - 1, face = "bold", hjust = 0.5))
}
nlp_max <- pairwise %>% mutate(nlp = if_else(is.finite(.data$p_adj_BH), -log10(pmax(.data$p_adj_BH, .Machine$double.xmin)), NA_real_)) %>%
  summarise(m = max(.data$nlp, na.rm = TRUE)) %>% pull(.data$m)
fill_limits <- c(0, nlp_max)
h_sp <- make_pair_plot("spearman", "Spearman", fill_limits); h_f1 <- make_pair_plot("max_f1", "Max F1", fill_limits)
h_gain <- make_pair_plot("pauc_gain", "pAUC gain", fill_limits); h_loss <- make_pair_plot("pauc_loss", "pAUC loss", fill_limits)
heat_2x2 <- ((h_sp + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 8))) | h_f1) / (h_gain | h_loss) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.key.height = unit(0.20, "cm"), legend.key.width = unit(0.5, "cm"))

# ---- (d) per-tool accuracy vs tumor purity (Spearman) ---------------
pur_raw2regid <- c(
  "InferCNV_HMM17_TumorFreeNormal_CNV"="infercnv_cnv","InferCNV_HMM17_TumorFreeNormal_Expr"="infercnv_expr",
  "InferCNVPy_TumorFreeNormal"="infercnvpy","Numbat_TumorFreeNormal_CNV"="numbat_cnv","Numbat_TumorFreeNormal_Expr"="numbat_expr",
  "CONICSmat_TumorFreeNormal"="conicsmat","CopyKAT_TumorFreeNormal"="copykat","SCEVAN_TumorFreeNormal"="scevan",
  "XClone_TumorFreeNormal"="xclone","CalicoST"="calicost","CopyVAE"="copyvae","STARCH"="starch")
pur_raw2lab <- c(
  "InferCNV_HMM17_TumorFreeNormal_CNV"="InferCNV-CNV","InferCNV_HMM17_TumorFreeNormal_Expr"="InferCNV-Expr",
  "InferCNVPy_TumorFreeNormal"="InferCNVpy","Numbat_TumorFreeNormal_CNV"="Numbat-CNV","Numbat_TumorFreeNormal_Expr"="Numbat-Expr",
  "CONICSmat_TumorFreeNormal"="CONICSmat","CopyKAT_TumorFreeNormal"="CopyKAT","SCEVAN_TumorFreeNormal"="SCEVAN",
  "XClone_TumorFreeNormal"="XClone","CalicoST"="CalicoST","CopyVAE"="CopyVAE","STARCH"="STARCH")
pur_col <- setNames(t1$color, t1$id); PUR_HI <- c("InferCNV-CNV", "Numbat-CNV", "Numbat-Expr", "CalicoST"); PUR_GREY <- "#bdbdbd"
pur_dat <- read_csv(file.path(IN_SYN, "purity_bin_metrics_all_tools.csv")) %>%
  transmute(regid = unname(pur_raw2regid[.data$tool]), lab = unname(pur_raw2lab[.data$tool]), purity_bin = .data$purity_bin,
            purity_mid = as.numeric(sub("-.*$", "", .data$purity_bin)) + 0.05, spearman = .data$spearman) %>%
  mutate(is_hi = .data$lab %in% PUR_HI, col = ifelse(.data$is_hi, unname(pur_col[.data$regid]), PUR_GREY))
pur_med <- pur_dat %>% group_by(.data$purity_bin, .data$purity_mid) %>% summarise(spearman = median(.data$spearman), .groups = "drop") %>% mutate(lab = "median", col = "black")
PUR_BINLABS <- pur_dat %>% distinct(.data$purity_mid, .data$purity_bin) %>% arrange(.data$purity_mid)
pur_others <- pur_dat %>% filter(!.data$is_hi); pur_hi <- pur_dat %>% filter(.data$is_hi)
pur_ends <- bind_rows(pur_hi %>% group_by(.data$lab) %>% slice_max(.data$purity_mid, n = 1) %>% ungroup(), pur_med %>% slice_max(.data$purity_mid, n = 1))
pD_purity <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.05) +
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, color = "grey40") +
  geom_line(data = pur_others, aes(.data$purity_mid, .data$spearman, group = .data$lab, color = .data$col), linewidth = 0.4, alpha = 0.55) +
  geom_line(data = pur_med, aes(.data$purity_mid, .data$spearman, color = .data$col), linetype = "dashed", linewidth = 0.9) +
  geom_line(data = pur_hi, aes(.data$purity_mid, .data$spearman, group = .data$lab, color = .data$col), linewidth = 0.9, alpha = 0.6) +
  geom_point(data = pur_hi, aes(.data$purity_mid, .data$spearman, color = .data$col), size = 1.0, alpha = 0.6, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = pur_ends, aes(.data$purity_mid, .data$spearman, label = .data$lab, color = .data$col),
                           family = FONT, hjust = 0, direction = "y", nudge_x = 0.03, size = 2.0, segment.size = 0.2, min.segment.length = 0, xlim = c(0.99, NA), show.legend = FALSE) +
  scale_color_identity(guide = "legend", name = NULL, breaks = c("black", PUR_GREY), labels = c("median (all 12 tools)", "other tools (n=8)")) +
  scale_x_continuous(breaks = PUR_BINLABS$purity_mid, labels = PUR_BINLABS$purity_bin, expand = expansion(mult = c(0.03, 0.34))) +
  coord_cartesian(ylim = c(0.0, 0.8)) +
  labs(title = "Per-tool accuracy vs tumor purity (Spearman)", x = "Tumor proportion (bin)", y = "Spearman") +
  theme_pub(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, margin = margin(b = 1)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 2), panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
        legend.position = "right", legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = BASE - 2))
write_tab(pur_dat %>% transmute(tool = lab, regid, purity_bin, purity_mid, spearman, highlighted = is_hi), "figS_task1A_panel_d_purity_curve_spearman.csv")

# ---- compose --------------------------------------------------------
design <- "
AAAAAAAA
BBBBB#CC
BBBBB#CC
DDDDD#CC
"
final <- wrap_plots(
  A = pA_type + labs(tag = "a"),
  B = wrap_elements(full = heat_2x2),
  C = pB_scope + labs(tag = "c"),
  D = pD_purity + labs(tag = "d"),
  design = design, widths = c(1.25, 1.25, 1.25, 1.25, 1.25, 0.15, 1.25, 1.25), heights = c(1.5, 1.75, 1.75, 1.35)
) & theme(plot.tag = element_text(face = "bold", size = 8), text = element_text(family = FONT))

W <- 7.09; H <- 8.4
ggsave(file.path(OUT_FIG, "figS_task1A_syn_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "figS_task1A_syn_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")

# ---- per-panel export -----------------------------------------------
# Panel b exported large (7.09 in square): the 2x2 pairwise matrices have small
# cells, so a bigger canvas grows the cells and keeps the significance stars
# legible when the panel is viewed on its own (composed heat_2x2 is unchanged).
panels <- list(a = list(pA_type, 7.09, 2.6), b = list(heat_2x2, 7.09, 7.09), c = list(pB_scope, 2.6, 5.2), d = list(pD_purity, 4.6, 2.6))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("figS_task1A_panel_", nm))
  pp <- pl[[1]] & theme(text = element_text(family = FONT))
  ggsave(paste0(base, ".pdf"), pp, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pp, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] figS_task1A_syn_compose + panels a-d (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
