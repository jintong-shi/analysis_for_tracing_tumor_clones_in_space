#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig S_1B (Task 1, REAL all-sample accuracy heatmap) — rebuilt.
#
# Reproduces paper/figures/supp_task1_merge/figS_task1B_real_heatmap_compose.R
# using ONLY input_data/real. The original selected the best-case run per
# sample/tool from all_datasets.csv; BEST.csv is that frozen selection, so it is
# consumed directly here (22 samples x 13 tools x 4 metrics, missing cells "na").
#
#   a Best-case accuracy heatmap: x = sample, y = tool, facet rows = 4 metrics.
#
# Sources only utils/figure_setup.R. Outputs: fig (pdf/png) in supp_figures/,
# full plotted grid in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(readr); library(purrr); library(tibble)
})

.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]))); normalizePath("figS_task1B_real_heatmap.R")
})()
SCRIPTS_DIR <- dirname(.this_file); TASK_DIR <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS   <- normalizePath(file.path(TASK_DIR, "..", "utils", "figure_setup.R"))
IN_REAL <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG <- file.path(TASK_DIR, "supp_figures"); OUT_TAB <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(UTILS)
source(file.path(.figure_setup_dir, "final_record_select.R"))

FONT <- pick_font()   # real "Arial" if installed, else bundled "Liberation Sans" (Arial metric clone), else "sans"
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

BASE <- 7
METRICS4 <- c("Spearman", "max_F1", "AUC_gain_trunc", "AUC_loss_trunc")
METRIC_LABELS <- c(Spearman = "Spearman", max_F1 = "Max F1",
                   AUC_gain_trunc = "pAUC gain", AUC_loss_trunc = "pAUC loss")
CNV_NEGATIVE <- c("skin_P4_V1", "skin_P4_V2")

.alias_lc <- tibble(key = tolower(TOOL_ALIAS_TABLE$alias), id = TOOL_ALIAS_TABLE$id) %>%
  distinct(key, .keep_all = TRUE)
norm_tool_id <- function(x) { id <- .alias_lc$id[match(tolower(as.character(x)), .alias_lc$key)]
  ifelse(is.na(id), as.character(x), id) }

best0 <- final_record_select_task1(read_csv(file.path(IN_REAL, "pseudobulk_vs_wes_metrics_all_datasets.csv"))) %>%
  mutate(sample = trimws(as.character(.data$sample)),
         tool_id = norm_tool_id(.data$tool),
         tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab),
         dataset_lab = dataset_label(as.character(.data$dataset))) %>%
  filter(!.data$sample %in% CNV_NEGATIVE)

sample_order_tbl <- best0 %>% distinct(.data$dataset, .data$dataset_lab, .data$sample) %>%
  arrange(.data$dataset_lab, .data$sample) %>% mutate(sample_x = .data$sample, x_num = row_number())
sample_order <- sample_order_tbl$sample_x
dataset_breaks <- sample_order_tbl %>% group_by(.data$dataset_lab) %>%
  summarise(last_x = max(.data$x_num), .groups = "drop") %>% arrange(.data$last_x) %>%
  filter(.data$last_x < length(sample_order)) %>% mutate(x = .data$last_x + 0.5)

present_meta <- TOOL_TABLE %>% filter(.data$id %in% unique(best0$tool_id)) %>% arrange(.data$order)
tool_order <- present_meta$label
metric_levels <- unname(METRIC_LABELS[METRICS4])

heat_values <- best0 %>%
  pivot_longer(all_of(METRICS4), names_to = "metric", values_to = "value") %>%
  mutate(metric = unname(METRIC_LABELS[.data$metric]),
         sample_x = as.character(.data$sample), tool_lab = as.character(.data$tool_lab)) %>%
  select("metric", "sample_x", "tool_lab", "value")

heat <- tidyr::expand_grid(metric = metric_levels, sample_x = sample_order, tool_lab = tool_order) %>%
  left_join(heat_values, by = c("metric", "sample_x", "tool_lab")) %>%
  mutate(metric = factor(.data$metric, levels = metric_levels),
         sample_x = factor(.data$sample_x, levels = sample_order),
         tool_lab = factor(.data$tool_lab, levels = rev(tool_order)),
         label = if_else(is.finite(.data$value), sprintf("%.2f", .data$value), "na"),
         label_color = case_when(!is.finite(.data$value) ~ "grey45",
                                 .data$value <= 0.05 ~ "white", .data$value >= 0.90 ~ "white",
                                 TRUE ~ "grey15"))
fill_low <- min(0, floor(min(heat$value, na.rm = TRUE) * 10) / 10)

p <- ggplot(heat, aes(.data$sample_x, .data$tool_lab, fill = .data$value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = .data$label, color = .data$label_color), size = 1.0, fontface = "bold") +
  geom_vline(data = dataset_breaks, aes(xintercept = .data$x), inherit.aes = FALSE, color = "grey25", linewidth = 0.25) +
  scale_x_discrete(labels = function(x) sub("^skin_", "SCC_", x)) +
  facet_grid(metric ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0.5,
                       limits = c(fill_low, 1), oob = scales::squish, na.value = "grey90", name = "Metric\nvalue") +
  scale_color_identity() +
  labs(title = "Real accuracy heatmap: best-case sample x tool values", x = "Sample (grouped by dataset)", y = NULL) +
  theme_pub(base_size = BASE) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, size = BASE - 2), axis.text.y = element_text(size = BASE - 2),
        strip.text.y = element_text(size = BASE - 1, face = "bold", angle = 0),
        legend.position = "right", legend.key.size = unit(0.28, "cm"),
        plot.title = element_text(size = BASE, face = "bold", hjust = 0))
final <- p + labs(tag = "a") & theme(plot.tag = element_text(face = "bold", size = 9), text = element_text(family = FONT))

W <- 7.09; H <- 8.2
ggsave(file.path(OUT_FIG, "figS_task1B_real_heatmap_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "figS_task1B_real_heatmap_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")
readr::write_csv(heat %>% arrange(.data$metric, .data$sample_x, .data$tool_lab),
                 file.path(OUT_TAB, "figS_task1B_real_heatmap_data.csv"))
cat(sprintf("[Done] figS_task1B_real_heatmap_compose (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
