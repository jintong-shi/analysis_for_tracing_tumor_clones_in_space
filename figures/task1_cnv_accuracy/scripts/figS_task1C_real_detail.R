#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Fig S_1C (Task 1, REAL detail panels) — rebuilt from scratch.
#
# Reproduces paper/figures/supp_task1_merge/figS_task1C_real_detail_compose.R
# using ONLY input_data/real. Panel letters (as composed):
#   a Ding FFPE pairwise paired-Wilcoxon significance (4 metrics)   [BEST]
#   b Ding FF pairwise paired-Wilcoxon significance (4 metrics)     [BEST]
#   c Ding per-class gain vs loss sensitivity by tool              [ding_gain_loss_*]
#   d SCC pt6 FF with BAM, Spearman                                 [BEST]
#   e Gene-level CNV heatmap, SCC_P6_V3 (WES + tools)               [pt6 WES + pseudobulk_gene]
#
# The tool rows of panel e read per-tool pseudobulk gene CSVs now stored under
# input_data/real/pseudobulk_gene/skin_P6_V3/ (see PROVENANCE_pseudobulk_gene.csv).
# Sources only utils/figure_setup.R. Outputs: composed fig + panels in
# supp_figures/; stats + plotted data in data_tables/.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(readr); library(purrr); library(tibble)
})

.this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]))); normalizePath("figS_task1C_real_detail.R")
})()
SCRIPTS_DIR <- dirname(.this_file); TASK_DIR <- normalizePath(file.path(SCRIPTS_DIR, ".."))
UTILS   <- normalizePath(file.path(TASK_DIR, "..", "utils", "figure_setup.R"))
IN_REAL <- file.path(TASK_DIR, "input_data", "real")
OUT_FIG <- file.path(TASK_DIR, "supp_figures"); OUT_TAB <- file.path(TASK_DIR, "data_tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
source(UTILS)
write_tab <- function(df, name) readr::write_csv(df, file.path(OUT_TAB, name))

FONT <- pick_font()   # real "Arial" if installed, else bundled "Liberation Sans" (Arial metric clone), else "sans"
update_geom_defaults("text",  list(family = FONT))
update_geom_defaults("label", list(family = FONT))

BASE <- 7
METRICS2 <- c("Spearman")
METRICS4 <- c("Spearman", "max_F1", "AUC_gain_trunc", "AUC_loss_trunc")
METRIC_LABELS <- c(Spearman = "Spearman", max_F1 = "Max F1",
                   AUC_gain_trunc = "pAUC gain", AUC_loss_trunc = "pAUC loss")
CNV_NEGATIVE <- c("skin_P4_V1", "skin_P4_V2")
STATE_TOOL_IDS <- c("infercnv_cnv", "numbat_cnv")   # S_1C heatmap: starch is log2 here
CHR_ORDER <- paste0("chr", 1:22); CN_LIMS <- c(1, 3)

.alias_lc <- tibble(key = tolower(TOOL_ALIAS_TABLE$alias), id = TOOL_ALIAS_TABLE$id) %>%
  distinct(key, .keep_all = TRUE)
norm_tool_id <- function(x) { id <- .alias_lc$id[match(tolower(as.character(x)), .alias_lc$key)]
  ifelse(is.na(id), as.character(x), id) }
metric_labeller <- as_labeller(function(x) { o <- unname(METRIC_LABELS[as.character(x)]); ifelse(is.na(o), as.character(x), o) })

source(file.path(.figure_setup_dir, "final_record_select.R"))
ALL_CSV <- file.path(IN_REAL, "pseudobulk_vs_wes_metrics_all_datasets.csv")
OVERVIEW_CSV <- file.path(.figure_setup_dir, "sample_overview.csv")
ov <- read_csv(OVERVIEW_CSV)
ov_sample <- trimws(as.character(ov$sample))

# best0 = BEST.csv (frozen best-run-per-(sample,tool) table).
best0 <- final_record_select_task1(read_csv(ALL_CSV)) %>%
  mutate(sample = trimws(as.character(.data$sample)), tool_id = norm_tool_id(.data$tool),
         tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab),
         dataset_lab = dataset_label(as.character(.data$dataset)),
         sample_type = ifelse(tolower(trimws(.data$sample_type)) == "ffpe", "ffpe", "ff")) %>%
  filter(.data$sample %in% ov_sample, !.data$sample %in% CNV_NEGATIVE)

present_meta <- TOOL_TABLE %>% filter(.data$id %in% unique(best0$tool_id)) %>% arrange(.data$order)
TOOL_ORDER <- present_meta$label; TOOL_COLORS_D <- setNames(present_meta$color, present_meta$label)
DING_LABEL <- dataset_label("Ding_nature_wes")

ding_ff_spearman_rank_ids <- best0 %>%
  filter(.data$dataset_lab == DING_LABEL, .data$sample_type == "ff", is.finite(.data$Spearman)) %>%
  group_by(.data$tool_id) %>% summarise(m = median(.data$Spearman, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(.data$m)) %>% arrange(desc(.data$m)) %>% pull(.data$tool_id)

# =====================================================================
# (d) SCC pt6 FF + BAM  (Spearman)
# =====================================================================
skin_pt6 <- best0 %>% filter(.data$sample %in% c("skin_P6_V3", "skin_P6_V4")) %>%
  pivot_longer(all_of(METRICS2), names_to = "metric", values_to = "value") %>%
  filter(is.finite(.data$value)) %>%
  mutate(metric = factor(.data$metric, levels = METRICS2), tool_lab = factor(.data$tool_lab, levels = TOOL_ORDER))
sample_pal <- setNames(c("#1F78B4", "#E31A1C"), sort(unique(skin_pt6$sample)))
p_pt6 <- ggplot(skin_pt6, aes(.data$tool_lab, .data$value)) +
  geom_boxplot(aes(fill = .data$tool_lab), width = 0.6, outlier.shape = NA, alpha = 0.55, color = "grey25", linewidth = 0.35) +
  geom_point(aes(color = .data$sample), position = position_jitter(width = 0.24, height = 0, seed = 1), size = 0.9, alpha = 0.9) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y", labeller = metric_labeller) +
  scale_fill_manual(values = TOOL_COLORS_D, drop = FALSE, guide = "none") +
  scale_color_manual(values = sample_pal, name = "Section", labels = function(x) sub("^skin_", "SCC_", x)) +
  labs(title = "SCC pt6 FF with BAM", x = NULL, y = "Metric value") +
  theme_pub(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE - 2),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3), legend.position = "top",
        legend.key.size = unit(0.24, "cm"), strip.text = element_text(size = BASE - 1, face = "bold"))
write_tab(skin_pt6 %>% transmute(sample, tool_lab, metric, value), "figS_task1C_panel_d_scc_pt6.csv")

# =====================================================================
# (e) SCC_P6_V3 gene-level heatmap
# =====================================================================
log2ratio_to_cn <- function(x) 2 * (2 ^ x); state_to_cn <- function(x) 2 + x
normalize_chr <- function(x) { x <- as.character(x); x <- sub("^CHR", "chr", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE); paste0("chr", x) }

gene_tbl <- read.delim(file.path(IN_REAL, "hgTables_hg38_gencode.txt"), sep = "\t",
                       header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_tbl <- gene_tbl[, c("name2", "chrom", "cdsStart", "cdsEnd")]
gene_tbl$mid <- (gene_tbl$cdsStart + gene_tbl$cdsEnd) %/% 2L
gene_tbl <- gene_tbl[gene_tbl$chrom %in% CHR_ORDER, ] %>%
  group_by(.data$name2) %>% slice_max(.data$cdsEnd - .data$cdsStart, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(chrom = factor(.data$chrom, levels = CHR_ORDER)) %>% arrange(.data$chrom, .data$mid) %>% mutate(gene_order = row_number())
chr_boundaries <- gene_tbl %>% group_by(.data$chrom) %>%
  summarise(first = min(.data$gene_order), last = max(.data$gene_order), .groups = "drop") %>%
  mutate(mid = (.data$first + .data$last) / 2)
chr_vlines <- chr_boundaries %>% mutate(end = .data$last + 0.5) %>% filter(.data$chrom != tail(CHR_ORDER, 1)) %>% pull(.data$end)

PSEUDO_DIR <- file.path(IN_REAL, "pseudobulk_gene", "skin_P6_V3")
TRUTH_WES_ID <- "P6_T"
load_pseudo <- function(path, tid) {
  d <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- "gene"; val_col <- setdiff(names(d), "gene")[1]
  raw_v <- suppressWarnings(as.numeric(d[[val_col]]))
  cn <- if (tid %in% STATE_TOOL_IDS) state_to_cn(raw_v) else log2ratio_to_cn(raw_v)
  data.frame(gene = as.character(d$gene), copy_number = cn, tool_id = tid, stringsAsFactors = FALSE)
}
tool_files <- list.files(PSEUDO_DIR, pattern = "\\.csv$", full.names = TRUE)
if (!length(tool_files)) stop("Panel e: no tool pseudobulk files under ", PSEUDO_DIR)
tool_long <- purrr::map_dfr(tool_files, function(p) load_pseudo(p, sub("\\.csv$", "", basename(p))))

seg <- read.csv(file.path(IN_REAL, "pt6_wes_segments_results.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
keep <- Reduce(`|`, lapply(intersect(c("sample", "wes_sample_id"), names(seg)),
                           function(cl) as.character(seg[[cl]]) == TRUTH_WES_ID))
seg <- seg[keep, , drop = FALSE]; seg$CONTIG <- normalize_chr(seg$CONTIG); seg <- seg[seg$CONTIG %in% CHR_ORDER, ]
if (!nrow(seg)) stop("No WES segments for skin_P6_V3 / ", TRUTH_WES_ID)
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
wes_long <- data.frame(gene = wes$name2, copy_number = log2ratio_to_cn(wes$log2ratio), tool_id = "WES", stringsAsFactors = FALSE)
all_long <- bind_rows(wes_long, tool_long) %>%
  inner_join(gene_tbl %>% select(.data$name2, .data$gene_order, .data$chrom), by = c("gene" = "name2"))
tool_labels <- setNames(TOOL_TABLE$label, TOOL_TABLE$id)
tool_ids_used <- TOOL_TABLE %>% filter(.data$id %in% setdiff(unique(all_long$tool_id), "WES")) %>% arrange(.data$order) %>% pull(.data$id)
ranked_tool_ids <- c(ding_ff_spearman_rank_ids[ding_ff_spearman_rank_ids %in% tool_ids_used],
                     setdiff(tool_ids_used, ding_ff_spearman_rank_ids))
row_levels <- c("WES", ranked_tool_ids)
row_labels <- setNames(c("WES (truth)", unname(tool_labels[ranked_tool_ids])), row_levels)
heat_skin <- all_long %>% mutate(row_label = unname(row_labels[.data$tool_id]),
                                 row_id = factor(.data$row_label, levels = rev(row_labels)))
p_heat <- ggplot(heat_skin, aes(.data$gene_order, .data$row_id, fill = .data$copy_number)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = chr_vlines, color = "grey55", linewidth = 0.18) +
  scale_fill_gradient2(low = "#2B5CB8", mid = "white", high = "#C9342C", midpoint = 2, limits = CN_LIMS,
                       breaks = c(1, 2, 3), labels = c("Loss", "Neutral", "Gain"), oob = scales::squish, na.value = "grey92", name = "CNV\nstate") +
  scale_x_continuous(breaks = chr_boundaries$mid, labels = sub("^chr", "", as.character(chr_boundaries$chrom)), expand = c(0, 0)) +
  labs(title = "Gene-level CNV heatmap: SCC_P6_V3", x = "Chromosome", y = NULL) +
  theme_pub(base_size = BASE) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = BASE - 1), axis.text.y = element_text(size = BASE - 1),
        legend.position = "right", legend.key.size = unit(0.28, "cm"))
write_tab(heat_skin %>% transmute(row_label, tool_id, gene, chrom, gene_order, copy_number),
          "figS_task1C_panel_e_heatmap_SCC_P6_V3.csv")

# =====================================================================
# (a,b) Ding pairwise significance heatmaps (BEST.csv)
# =====================================================================
DING_METRIC_LABELS <- c(Spearman = "Spearman", max_F1 = "Max F1", AUC_gain_trunc = "pAUC gain", AUC_loss_trunc = "pAUC loss")
best_ding <- final_record_select_task1(read_csv(ALL_CSV)) %>% filter(.data$dataset == "Ding_nature_wes") %>%
  mutate(tool_id = norm_tool_id(.data$tool), tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab),
         sample_type = ifelse(tolower(trimws(.data$sample_type)) == "ffpe", "ffpe", "ff"))
ding_tool_meta <- TOOL_TABLE %>% filter(.data$id %in% unique(best_ding$tool_id)) %>% arrange(.data$order)
DING_TOOL_ORDER <- ding_tool_meta$label
sig_label <- function(p) dplyr::case_when(is.na(p) ~ "", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "")
star_pw   <- function(p) dplyr::case_when(is.na(p) ~ "NA", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "ns")
real_pairwise_for_metric <- function(d, tools) {
  tools <- tools[tools %in% unique(d$tool_lab)]; if (length(tools) < 2) return(tibble())
  pairs <- utils::combn(tools, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    wide <- d %>% filter(.data$tool_lab %in% pair) %>% select(sample, tool_lab, value) %>% pivot_wider(names_from = tool_lab, values_from = value)
    x <- wide[[pair[1]]]; y <- wide[[pair[2]]]; ok <- is.finite(x) & is.finite(y)
    p_value <- if (sum(ok) < 2) NA_real_ else tryCatch(wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
    tibble(group1 = pair[1], group2 = pair[2], n_pairs = sum(ok), p_value = p_value)
  })
}
ding_pw <- purrr::map_dfr(c("ffpe", "ff"), function(st) {
  tools_st <- DING_TOOL_ORDER[DING_TOOL_ORDER %in% unique(best_ding$tool_lab[best_ding$sample_type == st])]
  purrr::map_dfr(names(DING_METRIC_LABELS), function(mt) {
    d <- best_ding %>% filter(.data$sample_type == st) %>% transmute(sample = .data$sample, tool_lab = .data$tool_lab, value = .data[[mt]]) %>% filter(is.finite(.data$value))
    res <- real_pairwise_for_metric(d, tools_st)
    if (nrow(res)) { res$sample_type <- st; res$metric <- mt }
    res
  })
}) %>% group_by(.data$sample_type, .data$metric) %>%
  mutate(p_adj_BH = p.adjust(.data$p_value, method = "BH"), signif = sig_label(.data$p_adj_BH)) %>% ungroup()
write_tab(ding_pw %>% arrange(.data$sample_type, .data$metric, .data$group1, .data$group2), "figS_task1C_ding_pairwise_wilcox.csv")

make_sym_real <- function(d, tools) {
  bind_rows(transmute(d, row_tool = .data$group1, col_tool = .data$group2, p_adj_BH = .data$p_adj_BH, signif = .data$signif),
            transmute(d, row_tool = .data$group2, col_tool = .data$group1, p_adj_BH = .data$p_adj_BH, signif = .data$signif),
            tibble(row_tool = tools, col_tool = tools, p_adj_BH = NA_real_, signif = "")) %>%
    mutate(row_tool = factor(.data$row_tool, levels = rev(tools)), col_tool = factor(.data$col_tool, levels = tools),
           neg_log10_p = if_else(is.finite(.data$p_adj_BH), -log10(pmax(.data$p_adj_BH, .Machine$double.xmin)), NA_real_))
}
ding_nlp_max <- ding_pw %>% mutate(nlp = if_else(is.finite(.data$p_adj_BH), -log10(pmax(.data$p_adj_BH, .Machine$double.xmin)), NA_real_)) %>%
  summarise(m = max(.data$nlp, na.rm = TRUE)) %>% pull(.data$m)
DING_FILL_LIMITS <- c(0, ding_nlp_max)
make_ding_panel <- function(stratum, title, fill_limits) {
  tools <- DING_TOOL_ORDER[DING_TOOL_ORDER %in% unique(best_ding$tool_lab[best_ding$sample_type == stratum])]
  heat_both <- purrr::map_dfr(names(DING_METRIC_LABELS), function(mt) {
    sub <- ding_pw %>% filter(.data$sample_type == stratum, .data$metric == mt)
    hs <- make_sym_real(sub, tools); hs$metric <- factor(mt, levels = names(DING_METRIC_LABELS)); hs
  })
  ggplot(heat_both, aes(.data$col_tool, .data$row_tool, fill = .data$neg_log10_p)) +
    geom_tile(color = "grey88", linewidth = 0.22) +
    geom_text(aes(label = .data$signif), size = 2.0, vjust = 0.75, color = "grey15", na.rm = TRUE) +
    facet_wrap(~ metric, nrow = 1, labeller = as_labeller(DING_METRIC_LABELS)) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = expression(-log[10](p[adj])), na.value = "grey95", limits = fill_limits, oob = scales::squish) +
    coord_equal() + labs(title = title, x = NULL, y = NULL) +
    theme_pub(base_size = BASE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 2), axis.text.y = element_text(size = BASE - 2),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(size = BASE - 1, face = "bold"), legend.key.size = unit(0.22, "cm"),
          legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
          plot.title = element_text(size = BASE, face = "bold", hjust = 0))
}
p_ffpe <- make_ding_panel("ffpe", paste0(DING_LABEL, " — FFPE"), DING_FILL_LIMITS)
p_ff   <- make_ding_panel("ff",   paste0(DING_LABEL, " — FF"),   DING_FILL_LIMITS)

# =====================================================================
# (c) Ding per-class gain vs loss sensitivity by tool
# =====================================================================
SENS_CSV <- file.path(IN_REAL, "ding_gain_loss_per_sample_sensitivity.csv")
SENS_STATS_CSV <- file.path(IN_REAL, "ding_gain_loss_tool_boxplot_stats.csv")
dir_tab <- figure_settings("direction", "synth_task1")
SENS_FILL <- c(Gain = dir_tab$color[dir_tab$id == "pauc_gain"], Loss = dir_tab$color[dir_tab$id == "pauc_loss"])
sens_stats <- read_csv(SENS_STATS_CSV) %>% filter(.data$n_paired_samples > 0) %>%
  mutate(tool_id = norm_tool_id(.data$tool), tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab), star = star_pw(.data$p_bh))
sens_tool_order <- TOOL_TABLE %>% filter(.data$id %in% unique(sens_stats$tool_id)) %>% arrange(.data$order) %>% pull(.data$label)
sens_paired <- read_csv(SENS_CSV) %>%
  mutate(tool_id = norm_tool_id(.data$tool), tool_lab = unname(TOOL_DISPLAY[.data$tool_id]),
         tool_lab = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab)) %>%
  filter(is.finite(.data$sensitivity_gain), is.finite(.data$sensitivity_loss), .data$tool_lab %in% sens_tool_order)
sens_long <- sens_paired %>%
  pivot_longer(c("sensitivity_gain", "sensitivity_loss"), names_to = "direction", values_to = "sensitivity") %>%
  mutate(direction = factor(ifelse(.data$direction == "sensitivity_gain", "Gain", "Loss"), levels = c("Gain", "Loss")),
         tool_lab = factor(.data$tool_lab, levels = sens_tool_order))
SENS_DODGE <- 0.72
sens_stars <- sens_stats %>% filter(.data$tool_lab %in% sens_tool_order) %>%
  mutate(tool_lab = factor(.data$tool_lab, levels = sens_tool_order), x = as.integer(.data$tool_lab))
p_sens <- ggplot(sens_long, aes(.data$tool_lab, .data$sensitivity, fill = .data$direction)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.75, color = "grey25", linewidth = 0.3, position = position_dodge(width = 0.72)) +
  geom_point(aes(color = .data$direction), show.legend = FALSE, size = 0.55, alpha = 0.45,
             position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.72, seed = 1)) +
  geom_segment(data = sens_stars, aes(x = .data$x - SENS_DODGE / 4, xend = .data$x + SENS_DODGE / 4, y = 1.02, yend = 1.02),
               inherit.aes = FALSE, color = "grey35", linewidth = 0.3) +
  geom_text(data = sens_stars, aes(x = .data$x, y = 1.055, label = .data$star), inherit.aes = FALSE, size = 2.0, color = "grey15", vjust = 0) +
  scale_fill_manual(values = SENS_FILL, name = "CNV class") + scale_color_manual(values = SENS_FILL, guide = "none") +
  scale_y_continuous(limits = c(-0.03, 1.13), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = paste0(DING_LABEL, " — gain vs loss"), x = NULL, y = "Sensitivity") +
  theme_pub(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE - 2), panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        legend.position = "top", legend.key.size = unit(0.24, "cm"), plot.title = element_text(size = BASE, face = "bold", hjust = 0))
write_tab(sens_long %>% transmute(sample, tool_lab, direction, sensitivity), "figS_task1C_panel_c_gain_loss_sensitivity.csv")

# ---- compose (panel letters a..e as in the original) ----------------
row_cd <- (p_sens + labs(tag = "c")) | (p_pt6 + labs(tag = "d"))
final <- (p_ffpe + labs(tag = "a")) / (p_ff + labs(tag = "b")) / row_cd / (p_heat + labs(tag = "e")) +
  plot_layout(heights = c(2.0, 2.0, 2.3, 1.5)) &
  theme(plot.tag = element_text(face = "bold", size = 9), text = element_text(family = FONT))

W <- 7.09; H <- 8.4
ggsave(file.path(OUT_FIG, "figS_task1C_real_detail_compose.pdf"), final, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "figS_task1C_real_detail_compose.png"), final, width = W, height = H, dpi = 300, bg = "white")

panels <- list(a = list(p_ffpe, 6.5, 2.4), b = list(p_ff, 6.5, 2.4), c = list(p_sens, 4.2, 3.0),
               d = list(p_pt6, 3.4, 3.0), e = list(p_heat, 7.09, 2.2))
for (nm in names(panels)) {
  pl <- panels[[nm]]; base <- file.path(OUT_FIG, paste0("figS_task1C_panel_", nm))
  pp <- pl[[1]] & theme(text = element_text(family = FONT))
  ggsave(paste0(base, ".pdf"), pp, width = pl[[2]], height = pl[[3]], limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pp, width = pl[[2]], height = pl[[3]], dpi = 300, bg = "white", limitsize = FALSE)
}
cat(sprintf("[Done] figS_task1C_real_detail_compose + panels a-e (%.2f x %.2f in) -> %s\n", W, H, OUT_FIG))
