#!/usr/bin/env Rscript
# =====================================================================
# Shared setup for the figures_final Task-3 figures (Fig 4 + Supp Task 3).
#
# Sourced by scripts/fig4_task3_main.R and scripts/figS_task3_supp.R AFTER
# utils/figure_setup.R (both live in figures_final/utils, so the panel scripts
# reference only utils/ code). It carries the pieces those two scripts share:
#   - task-3 tool / run / section orders and colours (from figure_settings.tsv)
#   - t3_load()        read the frozen best-case metrics + build sample groups
#   - main_slice()     best-run-per-tool filter (identity on the best-case csv)
#   - palette_table() / fill_gradientn()   TSV-driven continuous fills
#   - heat_theme() / theme_t3() / fmt_p()  shared theme + p-value formatter
#   - build_spatial()  the truth+tools spatial map (reused by 4d-f and Sc-Sd)
#
# All values are recomputed from input_data at run time; nothing is copied from
# a rendered figure. Depends only on symbols defined in figure_setup.R.
# =====================================================================

# best-case selection (shared): defines final_record_select_task3(), used by main_slice + t3_load
source(file.path(.figure_setup_dir, "final_record_select.R"))

FONT <- pick_font()
BASE <- 7
update_geom_defaults("text", list(family = FONT))
theme_t3 <- function() theme(text = element_text(family = FONT))

# ---- task-3 canonical orders / colours ------------------------------
TOOL_ORDER <- setting_order("real_task3")                 # tool ids, task-3 order
RUN_ORDER  <- figure_settings("run", "real_task3")$id
EXPR       <- c("copykat", "scevan", "copyvae", "tuscan") # expression-only tools
SEC_LV     <- c("FF", "FFPE")
section_type_colors <- setting_vector("section_type", "global", "color")
SEC_COL    <- setNames(unname(section_type_colors[c("fresh_frozen", "ffpe")]), SEC_LV)
sample_type_colors  <- setting_vector("sample_type", "global", "color")

# ---- shared helpers (verbatim from the frozen composes) -------------
main_slice <- final_record_select_task3   # shared best-run-per-tool filter (utils/final_record_select.R)
fmt_p <- function(p) ifelse(is.na(p), "p = NA", paste0("p = ", formatC(signif(p, 2), format = "g")))

palette_table <- function(scope) figure_settings("palette", scope) %>%
  mutate(value = suppressWarnings(as.numeric(.data$label)))
fill_gradientn <- function(scope, name, limits = c(0, 1)) {
  rows <- palette_table(scope) %>% arrange(.data$value)
  v <- rows$value; col <- rows$color
  if (min(v, na.rm = TRUE) > limits[1]) { v <- c(limits[1], v); col <- c(col[1], col) }
  if (max(v, na.rm = TRUE) < limits[2]) { v <- c(v, limits[2]); col <- c(col, col[length(col)]) }
  scale_fill_gradientn(colours = col, values = scales::rescale(v, from = limits),
                       limits = limits, oob = scales::squish, name = name, na.value = "grey90")
}
heat_theme <- function() theme_pub(BASE) + theme_t3() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.title = element_text(size = BASE - 1, face = "plain"),
        legend.text = element_text(size = BASE - 2), legend.key.size = unit(0.28, "cm"))

# ---- frozen best-case metrics + per-sample grouping -----------------
# section_type is deterministic: nature_prostate_ff = fresh-frozen, else FFPE.
# group is defined by the per-sample truth counts:
#   no tumor & no benign -> no_truth ; no tumor -> all_normal ;
#   no benign -> all_tumor ; else -> both.
t3_load <- function(input_dir) {
  df <- final_record_select_task3(read_csv(file.path(input_dir, "metrics_combined_relabelled.csv"))) %>%
    mutate(tool = factor(.data$tool, levels = TOOL_ORDER),
           run  = factor(.data$run,  levels = RUN_ORDER),
           section_type = ifelse(.data$dataset_id == "nature_prostate_ff", "fresh_frozen", "ffpe"))
  samp <- df %>% group_by(dataset_id, sample) %>%
    summarise(section_type = first(section_type), sample_type = first(sample_type),
              n_tumor  = max(n_truth_tumor_sample,  na.rm = TRUE),
              n_normal = max(n_truth_normal_sample, na.rm = TRUE),
              tumor_prop = { v <- tumor_proportion[run == "AllSpots_no-ref"]
                             if (length(v) && any(!is.na(v))) mean(v, na.rm = TRUE)
                             else mean(tumor_proportion, na.rm = TRUE) }, .groups = "drop") %>%
    mutate(group = case_when(n_tumor == 0 & n_normal == 0 ~ "no_truth",
                             n_tumor == 0 & n_normal >  0 ~ "all_normal",
                             n_tumor >  0 & n_normal == 0 ~ "all_tumor", TRUE ~ "both"),
           section_label = recode(section_type, ffpe = "FFPE", fresh_frozen = "FF",
                                  .default = section_type))
  dfg <- df %>% left_join(samp %>% select(sample, group, tumor_prop), by = "sample") %>%
    mutate(section_label = recode(section_type, ffpe = "FFPE", fresh_frozen = "FF",
                                  .default = section_type))
  list(df = df, samp = samp, dfg = dfg)
}

# ---- spatial map: ground truth + per-tool prediction ----------------
# Dual legend (Ground truth vs Prediction). Reused by Fig 4d-f and Supp Sc-Sd.
.tcol <- setting_vector("spatial_truth", "real_task3", "color")
.tlab <- setting_vector("spatial_truth", "real_task3", "label")
.pcol <- setting_vector("spatial_pred",  "real_task3", "color")
.plab <- setting_vector("spatial_pred",  "real_task3", "label")
ALL_COL <- c(pos = unname(.pcol["pos"]), neg = unname(.pcol["neg"]),
             other = unname(.tcol["other"]), unassigned = unname(.pcol["unassigned"]),
             na = unname(.pcol["na"]))
SP_DRAW <- c("other", "na", "unassigned", "neg", "pos")
.legend_grob <- function(ids, labs, cols, title) {
  d <- data.frame(x = 1, y = seq_along(ids), k = factor(ids, levels = ids))
  g <- ggplot(d, aes(x, y, color = k)) + geom_point(size = 1.6) +
    scale_color_manual(values = setNames(cols, ids), labels = labs, name = title, breaks = ids) +
    theme_void(base_size = BASE) +
    theme(legend.title = element_text(face = "plain", size = BASE - 1, family = FONT, hjust = 0),
          legend.text = element_text(size = BASE - 2, family = FONT), legend.key.size = unit(0.28, "cm"))
  cowplot::get_legend(g)
}
# returns list(plot=, data=) so the panel script can also save the ordered
# plotted points as a data table.
build_spatial <- function(input_dir, stem, title) {
  sp <- read_csv(file.path(input_dir, paste0("spatial_", stem, ".csv")))
  panel_levels <- sp %>% distinct(panel_order, panel) %>% arrange(panel_order) %>% pull(panel)
  sp <- sp %>% mutate(panel = factor(panel, levels = panel_levels),
                      class = factor(class, levels = names(ALL_COL))) %>%
    arrange(panel, factor(class, levels = SP_DRAW))
  xr <- range(sp$x); yr <- range(sp$y); cx <- mean(xr); cy <- mean(yr)
  half <- max(diff(xr), diff(yr)) / 2 * 1.04
  main <- ggplot(sp, aes(x, y, color = class)) +
    geom_point(size = 0.35, stroke = 0) + facet_wrap(~ panel, nrow = 1) +
    scale_color_manual(values = ALL_COL, guide = "none") +
    scale_x_continuous(limits = c(cx - half, cx + half), expand = c(0, 0)) +
    scale_y_continuous(limits = c(cy - half, cy + half), expand = c(0, 0), trans = "reverse") +
    coord_fixed(ratio = 1) + labs(title = title) + theme_void(base_size = BASE) + theme_t3() +
    theme(strip.text = element_text(face = "bold", size = BASE - 1, family = FONT, margin = margin(b = 1)),
          panel.spacing = unit(1.5, "pt"),
          plot.title = element_text(face = "bold", size = BASE, family = FONT, hjust = 0, margin = margin(l = 16, b = 2)))
  tl <- .legend_grob(c("pos", "neg", "other"), unname(.tlab[c("pos", "neg", "other")]),
                     ALL_COL[c("pos", "neg", "other")], "Ground truth")
  pl <- .legend_grob(c("pos", "neg", "unassigned", "na"), unname(.plab[c("pos", "neg", "unassigned", "na")]),
                     ALL_COL[c("pos", "neg", "unassigned", "na")], "Prediction")
  legends <- cowplot::plot_grid(tl, pl, ncol = 1, align = "v")
  plot <- wrap_plots(main, wrap_elements(full = legends), widths = c(length(panel_levels) * 1.15, 1.5))
  list(plot = plot, data = sp)
}
