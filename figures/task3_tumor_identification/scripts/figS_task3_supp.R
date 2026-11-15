#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig (Task 3) — composed supplementary figure, rebuilt.
#
# Reproduces paper/figures/supp_task3_merge/supp_task3_compose.R using ONLY the
# curated input_data of figures_final/task3_tumor_identification. Every plotted
# value is recomputed from the frozen best-case metrics / spatial csvs.
#
# Panels:
#   a Spearman(metric vs tumor proportion), both-class n=15    [metrics best-case]
#   b Pairwise paired-Wilcoxon (BH) on balanced accuracy, FF|FFPE [metrics best-case]
#   c Spatial P4_P13,  truth + 4 tools                         [spatial_P4_P13.csv]
#   d Spatial pt15_P2, truth + 4 tools                         [spatial_pt15_P2.csv]
#
# Scripts sourced: utils/figure_setup.R + utils/task3_setup.R (both in utils/).
# Outputs: composed supp fig + each panel (pdf/png) in supp_figures/, per-panel
# plotted data + stats in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(cowplot); library(purrr)
})

# ---- paths (this script's location -> project tree) -----------------
.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("figS_task3_supp.R")
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
source(file.path(UTILS_DIR, "task3_setup.R"))

write_tab <- function(df, name) readr::write_csv(df, file.path(OUT_TAB, name))

# ---- data + sample grouping (frozen best-case metrics) --------------
L    <- t3_load(IN_REAL)
df   <- L$df; samp <- L$samp; dfg <- L$dfg

# =====================================================================
# a  Spearman(metric vs tumor proportion) — all 20 samples
# =====================================================================
build_spearman <- function() {
  mets <- c("f1_score", "sensitivity", "specificity", "balanced_accuracy")
  met_lab <- c(f1_score = "F1", sensitivity = "Sens.", specificity = "Spec.", balanced_accuracy = "Bal. acc.")
  d4 <- dfg %>% filter(group == "both") %>% main_slice()   # both-class samples only (n=15)
  rho_tab <- d4 %>% dplyr::select(tool, tumor_prop, all_of(mets)) %>%
    tidyr::pivot_longer(all_of(mets), names_to = "metric", values_to = "val") %>%
    group_by(tool, metric) %>%
    summarise(rho = { ok <- !is.na(val) & !is.na(tumor_prop)
                      if (sum(ok) >= 4) suppressWarnings(cor(tumor_prop[ok], val[ok], method = "spearman")) else NA_real_ },
              p = { ok <- !is.na(val) & !is.na(tumor_prop)
                    if (sum(ok) >= 4) suppressWarnings(cor.test(tumor_prop[ok], val[ok], method = "spearman")$p.value) else NA_real_ },
              .groups = "drop") %>%
    mutate(tool_d = factor(unname(TOOL_DISPLAY[as.character(tool)]), levels = unname(TOOL_DISPLAY[rev(TOOL_ORDER)])),
           metric = factor(metric, levels = mets, labels = unname(met_lab[mets])),
           lab = ifelse(is.na(rho), "n/a", sprintf("%.2f", rho)),
           star = ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))),
           txt = ifelse(!is.na(rho) & abs(rho) > 0.65, "white", "black"))
  plot <- ggplot(rho_tab, aes(metric, tool_d, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = lab, color = txt), size = (BASE - 2.5) / .pt, nudge_y = 0.13, family = FONT) +
    geom_text(aes(label = star, color = txt), size = (BASE - 1.5) / .pt, nudge_y = -0.20, family = FONT) +
    scale_color_identity() + fill_gradientn("real_task3_spearman_rho", "Spearman\nrho", c(-1, 1)) +
    coord_equal() +
    labs(title = "Metric vs tumor proportion (Spearman)", x = NULL, y = NULL) + heat_theme() +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, margin = margin(l = 0, b = 2)),
          plot.title.position = "plot")
  list(plot = plot,
       data = rho_tab %>% transmute(tool = as.character(tool), metric = as.character(metric),
                                    spearman_rho = rho, p_value = p, star))
}

# =====================================================================
# b  Pairwise paired-Wilcoxon (BH) on balanced accuracy, FF | FFPE
# =====================================================================
build_pairwise <- function() {
  sig_star <- function(padj) ifelse(is.na(padj), "",
    ifelse(padj <= 0.001, "***", ifelse(padj <= 0.01, "**", ifelse(padj <= 0.05, "*", ""))))
  pairwise_ba <- function(dsub) {
    agg <- dsub %>% group_by(sample, tool) %>% summarise(BA = median(balanced_accuracy, na.rm = TRUE), .groups = "drop")
    tools <- TOOL_ORDER[TOOL_ORDER %in% as.character(unique(agg$tool))]
    if (length(tools) < 2) return(NULL)
    res <- purrr::map_dfr(utils::combn(tools, 2, simplify = FALSE), function(pr) {
      wide <- agg %>% filter(as.character(tool) %in% pr) %>% mutate(tool = as.character(tool)) %>%
        tidyr::pivot_wider(names_from = tool, values_from = BA)
      x <- wide[[pr[1]]]; y <- wide[[pr[2]]]; ok <- is.finite(x) & is.finite(y); n <- sum(ok)
      p <- if (n < 2) NA_real_ else tryCatch(wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE)$p.value,
                                             error = function(e) NA_real_)
      tibble(g1 = pr[1], g2 = pr[2], p_value = p)
    })
    res %>% mutate(p_adj_BH = p.adjust(p_value, "BH"))
  }
  make_sym <- function(res, tools) bind_rows(
    transmute(res, row = g1, col = g2, p_value, p_adj_BH),
    transmute(res, row = g2, col = g1, p_value, p_adj_BH),
    tibble(row = tools, col = tools, p_value = NA_real_, p_adj_BH = NA_real_)) %>%
    mutate(row_d = factor(unname(TOOL_DISPLAY[row]), levels = unname(TOOL_DISPLAY[rev(tools)])),
           col_d = factor(unname(TOOL_DISPLAY[col]), levels = unname(TOOL_DISPLAY[tools])),
           nlp = ifelse(is.finite(p_value), -log10(pmax(p_value, .Machine$double.xmin)), NA_real_),
           star = sig_star(p_adj_BH))
  d3 <- dfg %>% filter(group == "both") %>% main_slice() %>% filter(!is.na(balanced_accuracy)) %>%
    mutate(section_label = factor(section_label, levels = SEC_LV))
  res_ff <- pairwise_ba(d3 %>% filter(section_label == "FF"))
  res_pe <- pairwise_ba(d3 %>% filter(section_label == "FFPE"))
  nlp_max <- bind_rows(res_ff, res_pe) %>%
    mutate(nlp = ifelse(is.finite(p_value), -log10(pmax(p_value, .Machine$double.xmin)), NA_real_)) %>%
    summarise(m = max(nlp, na.rm = TRUE)) %>% pull(m)
  pw <- function(res, title) {
    tools <- TOOL_ORDER[TOOL_ORDER %in% unique(c(res$g1, res$g2))]
    ggplot(make_sym(res, tools), aes(col_d, row_d, fill = nlp)) +
      geom_tile(color = "grey88", linewidth = 0.3) +
      geom_text(aes(label = star), size = (BASE - 1.5) / .pt, vjust = 0.75, color = "grey15", na.rm = TRUE, family = FONT) +
      scale_fill_gradient(low = "white", high = "#B2182B", name = expression(-log[10](p)),
                          na.value = "grey95", limits = c(0, nlp_max), oob = scales::squish) +
      coord_equal() + labs(title = title, x = NULL, y = NULL) + heat_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = BASE, face = "bold", hjust = 0.5))
  }
  plot <- (pw(res_ff, "FF (n=5)") | pw(res_pe, "FFPE (n=10)")) + plot_layout(guides = "collect") +
    plot_annotation(title = "Pairwise paired-Wilcoxon (BH) on balanced accuracy",
      theme = theme(plot.title = element_text(face = "bold", size = BASE + 1, family = FONT,
                                              margin = margin(l = 16, b = 2))))
  tab <- bind_rows(res_ff %>% mutate(section = "FF"), res_pe %>% mutate(section = "FFPE")) %>%
    transmute(section, tool1 = g1, tool2 = g2, p_value, p_adj_BH, star = sig_star(p_adj_BH))
  list(plot = plot, stats = tab)
}

# =====================================================================
# BUILD panels + write data tables
# =====================================================================
r_a <- build_spearman()
r_b <- build_pairwise()
r_c <- build_spatial(IN_REAL, "P4_P13",  paste0(dataset_label("prostate_eva_ffpe"),       ", P4_P13"))
r_d <- build_spatial(IN_REAL, "pt15_P2", paste0(dataset_label("metastatic_prostate_ffpe"), ", pt15_P2"))

write_tab(r_a$data,  "figS_task3_panel_a_spearman_vs_tumorprop.csv")
write_tab(r_b$stats, "figS_task3_panel_b_pairwise_wilcox.csv")
write_tab(r_c$data %>% transmute(sample, panel, class = as.character(class), x, y),
          "figS_task3_panel_c_spatial_P4_P13.csv")
write_tab(r_d$data %>% transmute(sample, panel, class = as.character(class), x, y),
          "figS_task3_panel_d_spatial_pt15_P2.csv")

# =====================================================================
# COMPOSE + per-panel export
# =====================================================================
tp <- function(p) wrap_elements(full = p) + theme(plot.margin = margin(t = 4, r = 4, b = 2, l = 2))
top <- (tp(r_a$plot) | tp(r_b$plot)) + plot_layout(widths = c(1, 1.7))
final <- top / tp(r_c$plot) / tp(r_d$plot) +
  plot_layout(heights = c(2.6, 1.55, 1.55)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8, family = FONT),
        plot.tag.position = c(0, 1))

W <- 7.2; H <- 8.4
cairo_pdf(file.path(OUT_FIG, "supp_task3_compose.pdf"), width = W, height = H); print(final); dev.off()
png(file.path(OUT_FIG, "supp_task3_compose.png"), width = W, height = H, units = "in",
    res = 300, bg = "white", type = "cairo"); print(final); dev.off()

panels <- list(a = list(r_a$plot, 3.0, 3.2), b = list(r_b$plot, 4.4, 3.2),
               c = list(r_c$plot, 7.2, 1.9), d = list(r_d$plot, 7.2, 1.9))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("supp_task3_panel_", nm))
  cairo_pdf(paste0(base, ".pdf"), width = pl[[2]], height = pl[[3]]); print(pl[[1]]); dev.off()
  png(paste0(base, ".png"), width = pl[[2]], height = pl[[3]], units = "in",
      res = 300, bg = "white", type = "cairo"); print(pl[[1]]); dev.off()
}
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
cat(sprintf("[Done] supp_task3_compose + panels a-d (%.1f x %.1f in) -> %s\n", W, H, OUT_FIG))
