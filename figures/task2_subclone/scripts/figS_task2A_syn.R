#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig (Task 2, SYNTHETIC) — rebuilt from scratch using ONLY
# task2_subclone/input_data/synth.
#
# Panels:
#   a subclone-count overview heatmap (predicted - ground truth)  [subclone_ari_summary.csv]
#   b subclone-ARI pairwise Wilcoxon heatmap (multiclonal, GT>=2) [subclone_ari_summary.csv]
#   c n_clones sensitivity: subclone ARI vs requested n_clones     [subclone_ari_by_nclones.csv]
#
# Rule 6: the significance heatmap (b) colours AND stars are driven by the
# BH-adjusted p (p_adj_BH), not the raw p. The pairwise table is recomputed here
# from the summary (paired Wilcoxon per tool pair over shared reps, BH-adjusted).
# Sources ONLY the two utils scripts. Emits panels (pdf/png) in supp_figures/ and
# plotted data in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
  library(dplyr); library(stringr); library(scales)
})

.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("figS_task2A_syn.R")
})()
SCRIPTS_DIR <- dirname(.this_file)
TASK_DIR    <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS_DIR   <- normalizePath(file.path(TASK_DIR, "..", "utils"))
IN_SYN      <- file.path(TASK_DIR, "input_data", "synth")
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

BASE <- 7
sz   <- function(x) x * 0.78
pt   <- function(x) x * 0.72
theme_c <- function() {
  theme_task2_paper(base_size = BASE) +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0.5),
          plot.caption = element_text(size = BASE - 2, color = "gray40", hjust = 1))
}
sig_label <- function(p_adj) case_when(
  is.na(p_adj) ~ "", p_adj <= 0.001 ~ "***", p_adj <= 0.01 ~ "**", p_adj <= 0.05 ~ "*", TRUE ~ "")

SYN_CSV <- file.path(IN_SYN, "subclone_ari_summary.csv")
NCL_CSV <- file.path(IN_SYN, "subclone_ari_by_nclones.csv")
df <- load_subclone_summary(SYN_CSV)

# =====================================================================
# (a) subclone-count overview heatmap
# =====================================================================
rep_meta <- df %>%
  distinct(rep, n_gt_subclones) %>%
  mutate(rep_id = as.integer(str_extract(rep, "(?<=Rep_)\\d+")),
         rep_short = paste0("R", rep_id)) %>%
  arrange(n_gt_subclones, rep_id)
tools_long <- df %>%
  select(rep, tool, n_pred_subclones, status, n_gt_subclones) %>%
  mutate(pred_minus_gt = n_pred_subclones - n_gt_subclones,
         label = ifelse(status == "missing_prediction" | is.na(n_pred_subclones),
                        "—", as.character(n_pred_subclones)))
gt_long <- df %>%
  distinct(rep, n_gt_subclones) %>%
  mutate(tool = "Ground truth", n_pred_subclones = n_gt_subclones,
         pred_minus_gt = NA_real_, status = "ok", label = as.character(n_gt_subclones))
plot_df <- bind_rows(tools_long, gt_long) %>%
  left_join(rep_meta, by = c("rep", "n_gt_subclones")) %>%
  mutate(rep = factor(rep, levels = rep_meta$rep),
         rep_short = factor(rep_short, levels = rep_meta$rep_short),
         tool = factor(tool, levels = TOOL_ORDER_SYNTH_TASK2_OVERVIEW),
         text_white = !is.na(pred_minus_gt) & abs(pred_minus_gt) >= 10)
gt_yintercept <- length(TOOL_ORDER_SYNTH_TASK2_OVERVIEW) - 0.5
p_a <- ggplot(plot_df, aes(x = rep_short, y = tool, fill = pred_minus_gt)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(data = filter(plot_df, tool == "Ground truth"),
            fill = "gray85", color = "white", linewidth = 0.4) +
  geom_text(aes(label = label, color = text_white), size = sz(3.3), fontface = "bold") +
  geom_hline(yintercept = gt_yintercept, color = "black", linewidth = 0.5) +
  scale_fill_task2_pred_minus_gt(name = "Predicted\n− Ground truth") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  facet_grid(. ~ n_gt_subclones, scales = "free_x", space = "free_x",
             labeller = labeller(n_gt_subclones = function(x) paste0("GT=", x))) +
  labs(title = "Subclone count overview: predicted vs ground truth", x = NULL, y = NULL) +
  theme_task2_overview() +
  theme(plot.title = element_text(face = "bold", size = BASE + 1),
        axis.text.x = element_text(size = BASE - 2.5, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = BASE - 1),
        strip.text = element_text(size = BASE - 1, face = "bold"),
        legend.key.size = unit(0.30, "cm"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2))
write_tab(plot_df %>% filter(tool != "Ground truth") %>%
            transmute(rep, rep_short, tool, n_gt_subclones, n_pred_subclones, pred_minus_gt, status),
          "figS_task2A_panel_a_overview.csv")

# =====================================================================
# (b) subclone-ARI pairwise Wilcoxon heatmap (multiclonal, GT>=2)
#     Recomputed from the summary; BH-adjusted p drives colour AND stars.
# =====================================================================
multi <- df %>% filter(status == "ok", n_gt_subclones >= 2)
pw <- as.data.table(pairwise_paired_wilcox(
  data.frame(value = multi$ari, group = as.character(multi$tool), block = multi$rep)))
pw[, signif := sig_label(p_adj_BH)]
write_tab(pw[order(group1, group2)], "figS_task2A_panel_b_ari_pairwise_wilcox.csv")

present  <- unique(c(pw$group1, pw$group2))
tools_pw <- TOOL_ORDER_SYNTH_TASK2[TOOL_ORDER_SYNTH_TASK2 %in% present]
tools_pw <- c(tools_pw, setdiff(present, tools_pw))
heat <- rbind(
  pw[, .(row_tool = group1, col_tool = group2, p_adj_BH, signif)],
  pw[, .(row_tool = group2, col_tool = group1, p_adj_BH, signif)],
  data.table(row_tool = tools_pw, col_tool = tools_pw, p_adj_BH = NA_real_, signif = ""))
heat[, `:=`(row_tool = factor(row_tool, levels = rev(tools_pw)),
            col_tool = factor(col_tool, levels = tools_pw),
            neg_log10_padj = fifelse(is.finite(p_adj_BH),
                                     -log10(pmax(p_adj_BH, .Machine$double.xmin)), NA_real_))]
p_b <- ggplot(heat, aes(col_tool, row_tool, fill = neg_log10_padj)) +
  geom_tile(color = "grey88", linewidth = 0.25) +
  geom_text(aes(label = signif), size = 2.5, vjust = 0.76, color = "grey15", na.rm = TRUE) +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      name = expression(-log[10](p[adj])), na.value = "grey95") +
  coord_equal() +
  labs(title = "Subclone ARI pairwise Wilcoxon (multiclonal)", x = NULL, y = NULL) +
  theme_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = BASE - 1),
        panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right")

# =====================================================================
# (c) n_clones sensitivity: subclone ARI vs requested n_clones, per tool
#     Multiclone reps only; complete-case reps per tool (scored at EVERY
#     swept level) so boxplots + median + Friedman describe a balanced sample.
# =====================================================================
ncl <- as.data.frame(fread(NCL_CSV))
ncl$ari <- suppressWarnings(as.numeric(ncl$ari))
ncl <- ncl[ncl$n_gt_subclones > 1 & !is.na(ncl$ari), ]
ncl <- do.call(rbind, lapply(split(ncl, ncl$tool), function(sub) {
  nlev <- length(unique(sub$n_clones_param))
  keep <- names(which(table(sub$rep) == nlev))
  sub[sub$rep %in% keep, , drop = FALSE]
}))
ncl$tool <- factor(tool_label(as.character(ncl$tool)), levels = tool_label(c("STARCH", "CalicoST")))
ncl$n_clones_param <- factor(ncl$n_clones_param)
write_tab(ncl[, c("tool", "rep", "n_clones_param", "n_gt_subclones", "ari")],
          "figS_task2A_panel_c_nclones.csv")

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
ann <- do.call(rbind, lapply(split(ncl, ncl$tool), function(sub) {
  w  <- tapply(sub$ari, list(sub$rep, sub$n_clones_param), function(x) x[1])
  w  <- w[stats::complete.cases(w), , drop = FALSE]
  ft <- stats::friedman.test(w)
  data.frame(tool = sub$tool[1],
             n_clones_param = factor(tail(levels(ncl$n_clones_param), 1), levels = levels(ncl$n_clones_param)),
             ari = 1.0, statistic = unname(ft$statistic), df = unname(ft$parameter),
             p_value = ft$p.value, n_reps = nrow(w),
             label = sprintf("Friedman %s", fmt_p(ft$p.value)), stringsAsFactors = FALSE)
}))
write_tab(ann[, c("tool", "statistic", "df", "p_value", "n_reps")], "figS_task2A_panel_c_friedman.csv")
p_c <- ggplot(ncl, aes(n_clones_param, ari)) +
  geom_boxplot(outlier.shape = NA, fill = "#9ecae1", colour = "grey30", width = 0.6, alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, size = pt(1.5), alpha = 0.5, colour = "grey20") +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "#D55E00", linewidth = 0.7) +
  stat_summary(aes(group = 1), fun = median, geom = "point", colour = "#D55E00", size = pt(2.4)) +
  geom_text(data = ann, aes(label = label), hjust = 1, vjust = 1, size = sz(2.8)) +
  facet_wrap(~ tool, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(title = "n_clones sensitivity: subclone ARI vs requested n_clones",
       x = "n_clones parameter", y = "Subclone ARI") +
  theme_c() +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = BASE))

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
final <- (p_a + labs(tag = "a")) / (p_b + labs(tag = "b")) / (p_c + labs(tag = "c")) +
  plot_layout(heights = c(3.0, 2.4, 2.2)) &
  theme(plot.tag = element_text(face = "bold", size = 9), text = element_text(family = FONT))
W <- 7.09; H <- 8.0
ggsave(file.path(OUT_FIG, "figS_task2A_syn_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "figS_task2A_syn_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")

panels <- list(a = list(p_a, 7.09, 3.1), b = list(p_b, 3.6, 3.6), c = list(p_c, 4.6, 2.8))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("figS_task2A_panel_", nm))
  pp <- pl[[1]] & theme(text = element_text(family = FONT))
  ggsave(paste0(base, ".pdf"), pp, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pp, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] figS_task2A_syn_compose + panels a-c (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
cat("[Info] panel b pairwise (BH):\n"); print(pw[order(p_adj_BH), .(group1, group2, n_pairs, p_value, p_adj_BH, signif)])
cat("[Info] panel c Friedman:\n"); print(ann[, c("tool", "statistic", "df", "p_value", "n_reps")])
