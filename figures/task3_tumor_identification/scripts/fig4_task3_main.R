#!/usr/bin/env Rscript
# =====================================================================
# Fig 4 (Task 3: tumor-spot detection vs histology) — MAIN figure, rebuilt.
#
# Reproduces paper/figures/fig4_merge/fig4_compose.R using ONLY the curated
# input_data of figures_final/task3_tumor_identification. Every plotted value is
# recomputed from the frozen best-case metrics / spatial csvs; nothing is copied
# from a rendered figure.
#
# Panels:
#   a Sample overview + sample x tool detection score          [metrics best-case]
#   b Balanced accuracy by tool, FF | FFPE (both-class)        [metrics best-case]
#   c Negative controls: specificity by tool (all-normal)      [metrics best-case]
#   d Spatial V1_2  (all-normal),      truth + 5 tools         [spatial_V1_2.csv]
#   e Spatial H2_2  (normal-dominant), truth + 5 tools         [spatial_H2_2.csv]
#   f Spatial H1_5  (all-tumor),       truth + 5 tools         [spatial_H1_5.csv]
#
# Scripts sourced: utils/figure_setup.R + utils/task3_setup.R (both in utils/).
# Outputs: composed fig4 + each panel (pdf/png) in main_figures/, per-panel
# plotted data + stats in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(cowplot)
})

# ---- paths (this script's location -> project tree) -----------------
.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("fig4_task3_main.R")
})()
SCRIPTS_DIR <- dirname(.this_file)
TASK_DIR    <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS_DIR   <- normalizePath(file.path(TASK_DIR, "..", "utils"))
IN_REAL     <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG     <- file.path(TASK_DIR, "main_figures")
OUT_TAB     <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS_DIR, "figure_setup.R"))
source(file.path(UTILS_DIR, "task3_setup.R"))

write_tab <- function(df, name) readr::write_csv(df, file.path(OUT_TAB, name))

# ---- data + sample grouping (frozen best-case metrics) --------------
L    <- t3_load(IN_REAL)
df   <- L$df; samp <- L$samp; dfg <- L$dfg

# =====================================================================
# a  Sample overview + sample x tool detection score
# =====================================================================
build_overview <- function() {
  ord_all <- bind_rows(
    samp %>% filter(group == "all_normal") %>% arrange(tumor_prop),
    samp %>% filter(group == "both", section_label == "FF")   %>% arrange(tumor_prop),
    samp %>% filter(group == "both", section_label == "FFPE") %>% arrange(tumor_prop),
    samp %>% filter(group == "all_tumor") %>% arrange(tumor_prop))
  s_lev <- ord_all$sample
  GROUP_LV  <- c("all_normal", "both", "all_tumor")
  GROUP_COL <- setNames(unname(sample_type_colors[c("all_normal", "balanced", "all_tumor")]), GROUP_LV)
  METRIC_OF <- c(all_normal = "Specificity", both = "Balanced accuracy", all_tumor = "Sensitivity")
  ds_ids <- unique(as.character(ord_all$dataset_id))
  ds_col <- DATASET_COLORS[ds_ids]; ds_col[is.na(ds_col)] <- "grey60"; names(ds_col) <- ds_ids
  ds_lab <- dataset_label(ds_ids); names(ds_lab) <- ds_ids
  anno <- ord_all %>% mutate(sample = factor(sample, levels = s_lev),
                             group = factor(group, levels = GROUP_LV),
                             section_label = factor(section_label, levels = SEC_LV),
                             metric = METRIC_OF[as.character(group)])
  ba <- dfg %>% filter(group %in% GROUP_LV) %>% main_slice() %>%
    group_by(sample, tool, group) %>%
    summarise(spec = median(specificity, na.rm = TRUE), sens = median(sensitivity, na.rm = TRUE),
              bacc = median(balanced_accuracy, na.rm = TRUE), .groups = "drop") %>%
    mutate(BA = case_when(group == "all_normal" ~ spec, group == "all_tumor" ~ sens, TRUE ~ bacc),
           BA = ifelse(is.nan(BA), NA_real_, BA)) %>%
    dplyr::select(sample, tool, BA) %>%
    tidyr::complete(sample = s_lev, tool = TOOL_ORDER) %>%
    mutate(sample = factor(sample, levels = s_lev),
           tool_d = factor(unname(TOOL_DISPLAY[as.character(tool)]), levels = unname(TOOL_DISPLAY[rev(TOOL_ORDER)])),
           lab = ifelse(is.na(BA), "", sprintf("%.2f", BA)),
           txt = ifelse(!is.na(BA) & (BA < 0.4 | BA > 0.78), "white", "black"))
  strip_base <- function(p) p + labs(x = NULL, y = NULL) +
    scale_y_discrete(expand = expansion(add = 0)) + theme_pub(BASE) + theme_t3() +
    theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_text(size = BASE - 2),
          legend.position = "right", legend.title = element_text(size = BASE - 1, face = "plain"),
          legend.text = element_text(size = BASE - 2), legend.key.size = unit(0.28, "cm"),
          plot.margin = margin(0, 0, 0, 0))
  p_ch  <- strip_base(ggplot(anno, aes(sample, "Cohort", fill = as.character(dataset_id))) +
    geom_tile(color = "white", linewidth = 0.5) + scale_fill_manual(values = ds_col, labels = ds_lab, name = "Cohort"))
  p_tp  <- strip_base(ggplot(anno, aes(sample, "Tumor prop.", fill = tumor_prop)) +
    geom_tile(color = "white", linewidth = 0.5) + fill_gradientn("real_task3_tumor_proportion", "Tumor prop.", c(0, 1)))
  p_st  <- strip_base(ggplot(anno, aes(sample, "Sample type", fill = group)) +
    geom_tile(color = "white", linewidth = 0.5) + scale_fill_manual(values = GROUP_COL, name = "Sample type"))
  p_sec <- strip_base(ggplot(anno, aes(sample, "Section", fill = section_label)) +
    geom_tile(color = "white", linewidth = 0.5) + scale_fill_manual(values = SEC_COL, name = "Section"))
  p_met <- strip_base(ggplot(anno, aes(sample, "Metric", fill = metric)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(data = anno %>% group_by(metric, group) %>%
                summarise(sample = s_lev[round(mean(match(sample, s_lev)))], .groups = "drop") %>%
                mutate(hjust = ifelse(group == "all_tumor", 1, 0.5)),
              aes(sample, "Metric", label = metric, hjust = hjust), inherit.aes = FALSE,
              size = (BASE - 1) / .pt, family = FONT) +
    scale_fill_manual(values = c(Specificity = "#C7DCEF", `Balanced accuracy` = "#EDE6C8",
                                 Sensitivity = "#E7C9C9"), guide = "none"))
  p_main <- ggplot(ba, aes(sample, tool_d, fill = BA)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = lab, color = txt), size = (BASE - 2.5) / .pt) + scale_color_identity() +
    fill_gradientn("real_task3_balanced_accuracy", "Score\n(0-1)", c(0, 1)) +
    labs(x = NULL, y = NULL) + heat_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 2))
  plot <- p_ch / p_tp / p_st / p_sec / p_met / p_main +
    plot_layout(heights = c(1, 1, 1, 1, 1, 9), guides = "collect") +
    plot_annotation(title = "Sample overview & sample x tool detection score",
      theme = theme(plot.title = element_text(face = "bold", size = BASE + 1, family = FONT,
                                              margin = margin(l = 16, b = 2)))) &
    theme(legend.title = element_text(size = BASE - 1, face = "plain", family = FONT),
          legend.text = element_text(size = BASE - 2, family = FONT))
  list(plot = plot,
       score = ba %>% transmute(sample, tool = as.character(tool),
                                tool_label = as.character(tool_d), score = BA),
       meta  = anno %>% transmute(sample = as.character(sample), dataset_id = as.character(dataset_id),
                                  dataset_label = unname(ds_lab[as.character(dataset_id)]),
                                  section = as.character(section_label), group = as.character(group),
                                  tumor_prop, scoring_metric = metric))
}

# =====================================================================
# b  Balanced accuracy by tool, FF | FFPE (both-class samples)
# =====================================================================
build_ba_box <- function() {
  d2 <- dfg %>% filter(group == "both") %>% main_slice() %>% filter(!is.na(balanced_accuracy)) %>%
    mutate(tool = factor(tool, levels = TOOL_ORDER), section_label = factor(section_label, levels = SEC_LV),
           tool_d = factor(unname(TOOL_DISPLAY[as.character(tool)]), levels = unname(TOOL_DISPLAY[TOOL_ORDER])))
  kw <- d2 %>% group_by(section_label) %>%
    summarise(p = { t <- droplevels(tool)
                    if (nlevels(t) >= 2) suppressWarnings(kruskal.test(balanced_accuracy ~ t)$p.value) else NA_real_ },
              .groups = "drop") %>% mutate(lab = fmt_p(p))
  plot <- ggplot(d2, aes(tool_d, balanced_accuracy)) +
    geom_boxplot(aes(fill = tool), outlier.shape = NA, width = 0.62, alpha = 0.35, color = "black", linewidth = 0.3) +
    geom_jitter(aes(color = tool), width = 0.14, size = 0.9, alpha = 0.75, stroke = 0.2) +
    geom_text(data = kw, aes(x = Inf, y = 1.14, label = lab), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1, size = 5.5 / .pt, family = FONT) +
    facet_wrap(~ section_label, nrow = 1) +
    scale_color_manual(values = TOOL_COLORS) + scale_fill_manual(values = TOOL_COLORS) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25)) +
    labs(title = "Balanced accuracy by tool", x = NULL, y = "Balanced accuracy") +
    theme_pub(BASE) + theme_t3() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = BASE, face = "bold"),
          plot.title = element_text(face = "bold", size = BASE + 1, margin = margin(l = 16, b = 2)))
  list(plot = plot,
       data = d2 %>% transmute(sample, section = as.character(section_label),
                               tool = as.character(tool), balanced_accuracy),
       stats = kw %>% transmute(section = as.character(section_label),
                                test = "Kruskal-Wallis", p_value = p))
}

# =====================================================================
# c  Negative controls (all-normal): specificity by tool
# =====================================================================
build_spec_box <- function() {
  d11 <- bind_rows(
    dfg %>% filter(group == "all_normal") %>% main_slice() %>% filter(!is.na(specificity)) %>%
      transmute(sample, tool = as.character(tool), specificity, section_label),
    dfg %>% filter(group == "all_normal", section_label == "FF") %>% distinct(sample, section_label) %>%
      mutate(tool = "numbat", specificity = 1)) %>%
    mutate(tool = factor(tool, levels = TOOL_ORDER),
           tool_d = factor(unname(TOOL_DISPLAY[as.character(tool)]), levels = unname(TOOL_DISPLAY[TOOL_ORDER])),
           section_label = factor(section_label, levels = SEC_LV))
  nb_x <- which(levels(d11$tool_d) == unname(TOOL_DISPLAY["numbat"]))
  plot <- ggplot(d11, aes(tool_d, specificity)) +
    geom_boxplot(aes(fill = tool), outlier.shape = NA, width = 0.6, alpha = 0.25, color = "grey35", linewidth = 0.3) +
    geom_point(aes(color = section_label), position = position_jitter(width = 0.14, height = 0, seed = 1),
               size = 1.1, alpha = 0.9) +
    annotate("text", x = nb_x, y = 1.065, label = "no CNV called", size = 5 / .pt, color = "grey45", family = FONT) +
    scale_fill_manual(values = TOOL_COLORS, guide = "none") +
    scale_color_manual(values = SEC_COL, name = "Section", drop = FALSE) +
    scale_x_discrete(drop = FALSE) + scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    labs(title = "Negative controls: specificity by tool", x = NULL, y = "Specificity") +
    theme_pub(BASE) + theme_t3() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = BASE - 1, face = "plain"),
          legend.text = element_text(size = BASE - 2), legend.key.size = unit(0.28, "cm"),
          plot.title = element_text(face = "bold", size = BASE + 1, margin = margin(l = 16, b = 2)))
  list(plot = plot,
       data = d11 %>% transmute(sample, section = as.character(section_label),
                                tool = as.character(tool), specificity,
                                imputed = (tool == "numbat")))
}

# =====================================================================
# BUILD panels + write data tables
# =====================================================================
DL  <- dataset_label("nature_prostate_ff")
r_a <- build_overview()
r_b <- build_ba_box()
r_c <- build_spec_box()
r_d <- build_spatial(IN_REAL, "V1_2", paste0(DL, ", V1_2 (all-normal)"))
r_e <- build_spatial(IN_REAL, "H2_2", paste0(DL, ", H2_2 (normal-dominant)"))
r_f <- build_spatial(IN_REAL, "H1_5", paste0(DL, ", H1_5 (all-tumor)"))

write_tab(r_a$score, "fig4_panel_a_sample_tool_score.csv")
write_tab(r_a$meta,  "fig4_panel_a_sample_overview.csv")
write_tab(r_b$data,  "fig4_panel_b_balanced_accuracy.csv")
write_tab(r_b$stats, "fig4_panel_b_kruskal_wallis.csv")
write_tab(r_c$data,  "fig4_panel_c_negctrl_specificity.csv")
write_tab(r_d$data %>% transmute(sample, panel, class = as.character(class), x, y),
          "fig4_panel_d_spatial_V1_2.csv")
write_tab(r_e$data %>% transmute(sample, panel, class = as.character(class), x, y),
          "fig4_panel_e_spatial_H2_2.csv")
write_tab(r_f$data %>% transmute(sample, panel, class = as.character(class), x, y),
          "fig4_panel_f_spatial_H1_5.csv")

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
tp <- function(p) wrap_elements(full = p) + theme(plot.margin = margin(t = 4, r = 4, b = 2, l = 2))
row_bc <- tp(r_b$plot) | tp(r_c$plot)
final <- tp(r_a$plot) / row_bc / tp(r_d$plot) / tp(r_e$plot) / tp(r_f$plot) +
  plot_layout(heights = c(3.9, 2.3, 1.5, 1.5, 1.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8, family = FONT),
        plot.tag.position = c(0, 1))

W <- 7.2; H <- 10.2
cairo_pdf(file.path(OUT_FIG, "fig4_compose.pdf"), width = W, height = H); print(final); dev.off()
png(file.path(OUT_FIG, "fig4_compose.png"), width = W, height = H, units = "in",
    res = 300, bg = "white", type = "cairo"); print(final); dev.off()

panels <- list(a = list(r_a$plot, 7.2, 3.6), b = list(r_b$plot, 3.6, 2.6), c = list(r_c$plot, 3.6, 2.6),
               d = list(r_d$plot, 7.2, 1.7), e = list(r_e$plot, 7.2, 1.7), f = list(r_f$plot, 7.2, 1.7))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("fig4_panel_", nm))
  cairo_pdf(paste0(base, ".pdf"), width = pl[[2]], height = pl[[3]]); print(pl[[1]]); dev.off()
  png(paste0(base, ".png"), width = pl[[2]], height = pl[[3]], units = "in",
      res = 300, bg = "white", type = "cairo"); print(pl[[1]]); dev.off()
}
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
cat(sprintf("[Done] fig4_compose + panels a-f (%.1f x %.1f in) -> %s\n", W, H, OUT_FIG))
