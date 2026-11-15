#!/usr/bin/env Rscript
# =====================================================================
# Shared setup + statistics for the figures_final Task-2 (subclone) figures.
#
# This is a utils helper. Task-2 figure scripts source TWO utils files, in order:
#     source(".../utils/figure_setup.R")          # settings API, tool naming, theme_pub
#     source(".../utils/task2_subclone_setup.R")  # this file
# It ASSUMES figure_setup.R was sourced first (it uses figure_settings(),
# setting_order(), tool_label(), tool_id(), TOOL_DISPLAY, tool_color()). It adds
# only the Task-2-specific pieces that figure_setup.R does not carry:
#   - tool_color_vector()
#   - TOOL_ORDER_SYNTH_TASK2 / TOOL_COLORS_SYNTH_TASK2 / *_OVERVIEW
#   - theme_task2_paper() / theme_task2_overview()
#   - scale_fill_task2_pred_minus_gt()
#   - load_subclone_summary()
#   - omnibus stats: block_omnibus() / omnibus_label() / omnibus_row()
#     (Friedman for complete blocks, Skillings-Mack for incomplete)
#   - pairwise_paired_wilcox()  (paired Wilcoxon per tool pair, BH-adjusted)
#   - real FFPE|FF-BAM helpers: t2_load_real_bestcase(), t2_facet_median_order(),
#     t2_group_omni(), t2_group_pairwise()
#
# Ported verbatim (only the internal Skillings-Mack renamed to .sm_matrix so it
# does not clobber figure_setup.R's vector-arg skillings_mack) from the frozen
# paper/figures/{figure_setup.R, task2_subclone_stats.R}.
# =====================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Required R package is not installed: data.table", call. = FALSE)
  library(data.table)
})

# best-case selection (shared): defines final_record_select_task2() used by the loader below
source(file.path(.figure_setup_dir, "final_record_select.R"))

# ---- tool colour vector (figure_setup.R only ships tool_color()) -----
tool_color_vector <- function(x) setNames(tool_color(x), as.character(x))

# ---- Task-2 tool orders / colours (read from figure_settings.tsv) ----
TOOL_ORDER_SYNTH_TASK2 <- unname(TOOL_DISPLAY[setting_order("synth_task2_summary")])
TOOL_COLORS_SYNTH_TASK2 <- tool_color_vector(TOOL_ORDER_SYNTH_TASK2)
TOOL_ORDER_SYNTH_TASK2_OVERVIEW <-
  unname(TOOL_DISPLAY[setting_order("synth_task2_overview_ggplot")])

# ---- Task-2 themes (verbatim from paper/figures/figure_setup.R) ------
theme_task2_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 0.5, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "gray90", linetype = "dashed", linewidth = 0.3)
    )
}

theme_task2_overview <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.4, "lines"),
      strip.background = element_rect(fill = "gray92", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 9, color = "gray30"),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "gray40", size = 9.5)
    )
}

# ---- predicted - ground-truth diverging fill (from settings palette) --
scale_fill_task2_pred_minus_gt <- function(...) {
  pal <- figure_settings("palette", "synth_task2_pred_minus_gt")
  pal$value <- suppressWarnings(as.numeric(pal$label))
  scale_fill_gradient2(
    low  = pal$color[pal$id == "low"],
    mid  = pal$color[pal$id == "mid"],
    high = pal$color[pal$id == "high"],
    midpoint = pal$value[pal$id == "mid"],
    na.value = "gray85",
    limits = range(pal$value, na.rm = TRUE),
    breaks = c(-3, 0, 5, 10, 15, 20),
    ...
  )
}

# ---- synthetic subclone-summary loader -------------------------------
load_subclone_summary <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      tool = factor(.data$tool, levels = TOOL_ORDER_SYNTH_TASK2),
      is_monoclonal = .data$n_gt_subclones == 1
    )
}

# =====================================================================
# Omnibus statistics (Friedman / Skillings-Mack) across tools by block.
# `d` is a data.frame with numeric `value`, and `group` (tool) + `block`.
# =====================================================================

# Moore-Penrose pseudo-inverse via SVD (Sigma is singular, rank g-1).
.ginv_svd <- function(M, tol = 1e-9) {
  s <- svd(M)
  pos <- s$d > tol * max(s$d)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# Skillings-Mack statistic on a block x group matrix (NA = missing).
.sm_matrix <- function(mat) {
  obs  <- !is.na(mat)
  keep <- rowSums(obs) >= 2L
  mat  <- mat[keep, , drop = FALSE]
  obs  <- obs[keep, , drop = FALSE]
  g <- ncol(mat); b <- nrow(mat)
  A <- setNames(numeric(g), colnames(mat))
  Sigma <- matrix(0, g, g, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(b)) {
    present <- which(obs[j, ])
    k <- length(present)
    r <- rank(mat[j, present], ties.method = "average")
    a <- sqrt(12 / (k + 1)) * (r - (k + 1) / 2)
    A[present] <- A[present] + a
    Sigma[present, present] <- Sigma[present, present] - 1
    diag(Sigma)[present] <- diag(Sigma)[present] + k
  }
  stat <- as.numeric(t(A) %*% .ginv_svd(Sigma) %*% A)
  df   <- g - 1L
  list(method = "Skillings-Mack",
       statistic = stat, parameter = df,
       p.value = pchisq(stat, df, lower.tail = FALSE),
       n_blocks = b, n_groups = g)
}

# Omnibus test: Friedman (complete blocks) vs Skillings-Mack (incomplete).
block_omnibus <- function(d) {
  d <- as.data.frame(d)
  d <- d[!is.na(d$block) & !is.na(d$group), c("value", "group", "block")]
  mat <- tapply(d$value, list(d$block, d$group), function(x) mean(x, na.rm = TRUE))
  mat[is.nan(mat)] <- NA
  complete <- all(!is.na(mat))
  if (complete) {
    dd <- d
    dd$group <- factor(dd$group)
    dd$block <- factor(dd$block)
    ft <- friedman.test(value ~ group | block, data = dd)
    return(list(method = "Friedman",
                statistic = unname(ft$statistic),
                parameter = unname(ft$parameter),
                p.value = ft$p.value,
                n_blocks = nrow(mat), n_groups = ncol(mat)))
  }
  .sm_matrix(mat)
}

omnibus_label <- function(res)
  paste0(res$method, ", p = ", scales::pvalue(res$p.value, accuracy = 0.001))

omnibus_row <- function(res, metric)
  data.frame(metric = metric, method = res$method,
             statistic = res$statistic, df = res$parameter,
             p_value = res$p.value, n_blocks = res$n_blocks,
             n_groups = res$n_groups, stringsAsFactors = FALSE)

# Pairwise paired Wilcoxon signed-rank over shared blocks, BH-adjusted.
pairwise_paired_wilcox <- function(d, p_adjust = "BH") {
  d <- as.data.frame(d)
  d <- d[!is.na(d$block) & !is.na(d$group) & is.finite(d$value), ]
  groups <- sort(unique(as.character(d$group)))
  combs <- utils::combn(groups, 2)
  rows <- lapply(seq_len(ncol(combs)), function(c) {
    g1 <- combs[1, c]; g2 <- combs[2, c]
    w <- merge(d[d$group == g1, c("block", "value")],
               d[d$group == g2, c("block", "value")],
               by = "block", suffixes = c("_1", "_2"))
    n <- nrow(w)
    if (n >= 1L) {
      tt <- tryCatch(wilcox.test(w$value_1, w$value_2, paired = TRUE, exact = FALSE),
                     error = function(e) NULL)
      p    <- if (is.null(tt)) NA_real_ else tt$p.value
      stat <- if (is.null(tt)) NA_real_ else unname(tt$statistic)
      m1 <- median(w$value_1); m2 <- median(w$value_2)
    } else { p <- NA_real_; stat <- NA_real_; m1 <- NA_real_; m2 <- NA_real_ }
    data.frame(group1 = g1, group2 = g2, n_pairs = n,
               median1 = m1, median2 = m2,
               statistic = stat, p_value = p, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj_BH <- p.adjust(out$p_value, method = p_adjust)
  out[order(out$p_adj_BH), ]
}

# =====================================================================
# REAL FFPE | FF-BAM best-case table helpers. Input = the pre-filtered
# best-case table (one row per dataset/sample/tool, with a `group` column).
# =====================================================================
t2_load_real_bestcase <- function(csv,
                                  group_levels = c("FFPE", "FF-BAM"),
                                  ffpe_drop    = c("Numbat", "CalicoST")) {
  raw <- final_record_select_task2(data.table::fread(csv))
  raw[, group := factor(group, levels = group_levels)]
  raw[, tool  := tool_label(tool)]
  best <- raw[status == "ok" & is.finite(z_score) & is.finite(assortativity_obs)]
  nc   <- raw[status %in% c("ok", "single_subclone") & is.finite(n_subclones)]
  lv   <- intersect(TOOL_DISPLAY, union(unique(best$tool), unique(nc$tool)))
  best[, tool := factor(tool, levels = lv)]
  nc[,   tool := factor(tool, levels = lv)]
  key  <- best[, .(k = paste(dataset, sample, tool), assortativity_obs)]
  list(best = best, nc = nc,
       best_ab      = best[!(group == "FFPE" & tool %in% ffpe_drop)],
       tool_levels  = lv,
       tool_pal     = tool_color_vector(lv),
       group_levels = group_levels,
       assort = function(dataset, sample, tool) {
         v <- key[k == paste(dataset, sample, tool_label(tool)), assortativity_obs]
         if (length(v)) v[1L] else NA_real_
       })
}

# Per-facet x order by descending median of `valcol`; adds `xid` + label map.
t2_facet_median_order <- function(d, valcol) {
  m <- d[, .(median = median(get(valcol), na.rm = TRUE)),
         by = .(group, tool)][order(group, -median)]
  m[, xid := paste(group, tool, sep = "___")]
  d2 <- data.table::copy(d)
  d2[, xid := factor(paste(group, tool, sep = "___"), levels = m$xid)]
  list(d = d2, labels = setNames(as.character(m$tool), m$xid), med = m)
}

# Per-group omnibus label table (one row per group level present in `d_ab`).
t2_group_omni <- function(d_ab, valcol, group_levels = c("FFPE", "FF-BAM")) {
  data.table::rbindlist(lapply(split(d_ab, droplevels(d_ab$group)), function(d) {
    lab <- tryCatch(omnibus_label(block_omnibus(data.frame(
      value = d[[valcol]], group = as.character(d$tool),
      block = paste(d$dataset, d$sample, sep = "::")))), error = function(e) "")
    data.table::data.table(group = factor(as.character(d$group[1L]), group_levels), label = lab)
  }), fill = TRUE)
}

# Per-group pairwise paired-Wilcoxon table (adds a `platform` column).
t2_group_pairwise <- function(d, valcol) {
  data.table::rbindlist(lapply(split(d, droplevels(d$group)), function(g) {
    res <- data.table::as.data.table(pairwise_paired_wilcox(data.frame(
      value = g[[valcol]], group = as.character(g$tool),
      block = paste(g$dataset, g$sample, sep = "::"))))
    if (!nrow(res)) return(NULL)
    res[, platform := as.character(g$group[1L])]
    res
  }), fill = TRUE)
}
