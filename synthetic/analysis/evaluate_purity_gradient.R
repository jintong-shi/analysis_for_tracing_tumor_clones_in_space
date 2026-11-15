#!/usr/bin/env Rscript
# Purity-gradient evaluation against the purity-diluted ground truth in
# st-genome_profile.parquet (GT[gene,spot] = prop*clone_CN + (1-prop)*2).
# Truth labels use a CN band (default gain > 2.25, loss < 1.75) so only detectable
# dilutions count as events, and a per-spot Spearman is related to that spot's
# continuous purity alongside the bin-level pooled F1/AUC.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(arrow)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
})
EPS <- 1e-12
SPLIT_TOOL_NAMES <- c("InferCNV", "Numbat")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
make_purity_bin <- function(x, min_p = 0.2, max_p = 1.0, by = 0.1) {
  breaks <- seq(min_p, max_p, by = by)
  if (tail(breaks, 1) < max_p) breaks <- c(breaks, max_p)
  labels <- paste0(
    formatC(breaks[-length(breaks)], format = "f", digits = 1),
    "-",
    formatC(breaks[-1], format = "f", digits = 1)
  )
  as.character(cut(x, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE))
}

read_pred_matrix <- function(csv_file) {
  df <- fread(csv_file, data.table = FALSE)
  if (ncol(df) < 2) return(NULL)
  spot_ids <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  num_df <- as.data.frame(
    lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE
  )
  mat <- as.matrix(num_df)
  rownames(mat) <- spot_ids
  storage.mode(mat) <- "numeric"
  mat
}

state_matrix_to_numeric <- function(state_mat) {
  out <- matrix(0, nrow = nrow(state_mat), ncol = ncol(state_mat),
                dimnames = dimnames(state_mat))
  s <- tolower(as.character(as.matrix(state_mat)))
  out[s %in% c("gain", "amp", "bamp", "hamp", "amplified", "amplification")] <- 1
  out[s %in% c("loss", "del", "bdel", "deletion")] <- -1
  out[s %in% c("neutral", "neu", "loh", "cnloh", "cn_loh")] <- 0
  storage.mode(out) <- "numeric"
  out
}

read_prediction_matrix <- function(csv_file, matrix_type = c("log2", "state")) {
  matrix_type <- match.arg(matrix_type)
  if (matrix_type == "state") {
    df_in <- read.csv(csv_file, row.names = 1, check.names = FALSE)
    return(state_matrix_to_numeric(as.matrix(df_in)))
  }
  read_pred_matrix(csv_file)
}

get_tool_eval_specs <- function(tool, tool_out_dir) {
  split_tool <- tool %in% SPLIT_TOOL_NAMES ||
    grepl("^InferCNV_", tool) ||
    grepl("^Numbat_", tool)

  if (split_tool) {
    specs <- list(
      list(
        label = paste0(tool, "_CNV"),
        file = file.path(tool_out_dir, "processed_state_spots_x_genes.csv"),
        matrix_type = "state"
      ),
      list(
        label = paste0(tool, "_Expr"),
        file = file.path(tool_out_dir, "processed_log2_spots_x_genes.csv"),
        matrix_type = "log2"
      )
    )
    return(Filter(function(x) file.exists(x$file), specs))
  }

  if (tool == "XClone" || grepl("^XClone_", tool)) {
    state_file <- file.path(tool_out_dir, "processed_state_spots_x_genes.csv")
    if (!file.exists(state_file)) return(list())
    return(list(list(label = tool, file = state_file, matrix_type = "state")))
  }

  csv_file <- file.path(tool_out_dir, "processed_log2_spots_x_genes.csv")
  if (!file.exists(csv_file)) return(list())
  list(list(label = tool, file = csv_file, matrix_type = "log2"))
}

read_gt_matrix <- function(parquet_file) {
  gt_df <- read_parquet(parquet_file)
  if ("index" %in% names(gt_df)) gt_df <- tibble::column_to_rownames(gt_df, "index")
  gt_mat <- t(as.matrix(gt_df))
  storage.mode(gt_mat) <- "numeric"
  gt_mat
}

rbind_union_cols <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  all_cols <- union(colnames(a), colnames(b))

  a2 <- matrix(NA_real_, nrow = nrow(a), ncol = length(all_cols),
               dimnames = list(rownames(a), all_cols))
  b2 <- matrix(NA_real_, nrow = nrow(b), ncol = length(all_cols),
               dimnames = list(rownames(b), all_cols))

  a2[, colnames(a)] <- a
  b2[, colnames(b)] <- b
  rbind(a2, b2)
}

write_matrix_csv <- function(mat, out_file) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  df <- data.frame(spot_id = rownames(mat), mat, check.names = FALSE)
  fwrite(df, out_file)
}

has_global_cache <- function(base_out_dir, tool) {
  tool_dir <- file.path(base_out_dir, tool)
  file.exists(file.path(tool_dir, "global_pred_mat.csv.gz")) &&
    file.exists(file.path(tool_dir, "global_gt_mat.csv.gz")) &&
    file.exists(file.path(tool_dir, "global_metadata.csv"))
}

load_global_cache <- function(base_out_dir, tool) {
  tool_dir <- file.path(base_out_dir, tool)
  pred <- read_pred_matrix(file.path(tool_dir, "global_pred_mat.csv.gz"))
  gt <- read_pred_matrix(file.path(tool_dir, "global_gt_mat.csv.gz"))
  meta <- fread(file.path(tool_dir, "global_metadata.csv"), data.table = FALSE)

  if (!"spot_id" %in% names(meta) && ncol(meta) > 0) {
    names(meta)[1] <- "spot_id"
  }
  if (!"spot_id" %in% names(meta)) {
    stop(sprintf("[%s] global_metadata.csv has no 'spot_id' column.", tool))
  }

  list(pred = pred, gt = gt, meta = meta)
}

# ------------------------------------------------------------------------------
# Metrics (same as 05/06)
# ------------------------------------------------------------------------------
macro_f1_3class <- function(pred_class, true_class, classes = c(-1L, 0L, 1L)) {
  f1s <- numeric(length(classes))
  names(f1s) <- as.character(classes)

  for (k in seq_along(classes)) {
    cls <- classes[k]

    tp <- sum(pred_class == cls & true_class == cls, na.rm = TRUE)
    fp <- sum(pred_class == cls & true_class != cls, na.rm = TRUE)
    fn <- sum(pred_class != cls & true_class == cls, na.rm = TRUE)

    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

    f1 <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_
    if (is.na(f1)) f1 <- 0

    f1s[k] <- f1
  }

  mean(f1s)
}

evaluate_optimal_f1_score_log2ratio <- function(pred_score, true_class, baseline = 0, max_candidates = 200) {
  ok <- is.finite(pred_score) & !is.na(true_class)
  pred_score <- pred_score[ok]
  true_class <- as.integer(true_class[ok])

  if (length(pred_score) < 50) {
    return(c(max_f1=NA, gain_thr=NA, sens_gain=NA, prec_gain=NA,
             loss_thr=NA, sens_loss=NA, prec_loss=NA))
  }

  scores <- sort(unique(pred_score))

  loss_scores <- scores[scores <= baseline]
  if (length(loss_scores) > max_candidates) {
    loss_scores <- loss_scores[round(seq(1, length(loss_scores), length.out = max_candidates))]
  }
  if (length(loss_scores) == 0) loss_scores <- baseline
  loss_scores <- unique(c(loss_scores, baseline))

  gain_scores <- scores[scores >= baseline]
  if (length(gain_scores) > max_candidates) {
    gain_scores <- gain_scores[round(seq(1, length(gain_scores), length.out = max_candidates))]
  }
  if (length(gain_scores) == 0) gain_scores <- baseline
  gain_scores <- unique(c(baseline, gain_scores))

  ord <- order(pred_score)
  s <- pred_score[ord]
  y <- true_class[ord]

  r <- rle(s)
  u <- r$values
  ends <- cumsum(r$lengths)
  prev_ends <- c(0L, ends[-length(ends)])

  classes <- c(-1L, 0L, 1L)
  group_counts <- matrix(0L, nrow = length(u), ncol = 3)
  colnames(group_counts) <- c("loss", "neutral", "gain")

  for (k in seq_along(classes)) {
    cc <- cumsum(y == classes[k])
    prev_counts <- numeric(length(prev_ends))
    has_prev <- prev_ends > 0L
    prev_counts[has_prev] <- cc[prev_ends[has_prev]]
    group_counts[, k] <- cc[ends] - prev_counts
  }

  cum_counts <- rbind(c(0, 0, 0), apply(group_counts, 2, cumsum))
  colnames(cum_counts) <- c("loss", "neutral", "gain")

  total_cls <- cum_counts[nrow(cum_counts), ]
  names(total_cls) <- c("loss", "neutral", "gain")

  le_loss <- findInterval(loss_scores, u)
  eq_loss <- le_loss > 0L
  eq_loss[eq_loss] <- u[le_loss[eq_loss]] == loss_scores[eq_loss]
  lt_loss <- le_loss - as.integer(eq_loss)
  loss_counts <- cum_counts[lt_loss + 1L, , drop = FALSE]

  le_gain <- findInterval(gain_scores, u)
  gain_counts <- sweep(cum_counts[le_gain + 1L, , drop = FALSE], 2, total_cls, FUN = function(x, t) t - x)

  best_f1 <- -Inf
  best_loss_thr <- NA_real_
  best_gain_thr <- NA_real_

  for (j in seq_along(gain_scores)) {
    g_thr <- gain_scores[j]
    valid_i <- which(loss_scores <= g_thr)
    if (length(valid_i) == 0) next

    gc <- gain_counts[j, ]

    for (i in valid_i) {
      lc <- loss_counts[i, ]

      tp_l <- lc["loss"]
      fp_l <- sum(lc) - tp_l
      fn_l <- total_cls["loss"] - tp_l

      prec_l <- if ((tp_l + fp_l) > 0) tp_l / (tp_l + fp_l) else NA_real_
      rec_l  <- if ((tp_l + fn_l) > 0) tp_l / (tp_l + fn_l) else NA_real_
      f1_l <- if (!is.na(prec_l) && !is.na(rec_l) && (prec_l + rec_l) > 0) {
        2 * prec_l * rec_l / (prec_l + rec_l)
      } else {
        0
      }

      tp_g <- gc["gain"]
      fp_g <- sum(gc) - tp_g
      fn_g <- total_cls["gain"] - tp_g

      prec_g <- if ((tp_g + fp_g) > 0) tp_g / (tp_g + fp_g) else NA_real_
      rec_g  <- if ((tp_g + fn_g) > 0) tp_g / (tp_g + fn_g) else NA_real_
      f1_g <- if (!is.na(prec_g) && !is.na(rec_g) && (prec_g + rec_g) > 0) {
        2 * prec_g * rec_g / (prec_g + rec_g)
      } else {
        0
      }

      nc <- total_cls - lc - gc

      tp_n <- nc["neutral"]
      fp_n <- sum(nc) - tp_n
      fn_n <- total_cls["neutral"] - tp_n

      prec_n <- if ((tp_n + fp_n) > 0) tp_n / (tp_n + fp_n) else NA_real_
      rec_n  <- if ((tp_n + fn_n) > 0) tp_n / (tp_n + fn_n) else NA_real_
      f1_n <- if (!is.na(prec_n) && !is.na(rec_n) && (prec_n + rec_n) > 0) {
        2 * prec_n * rec_n / (prec_n + rec_n)
      } else {
        0
      }

      macro_f1 <- mean(c(f1_l, f1_n, f1_g))

      if (macro_f1 > best_f1) {
        best_f1 <- macro_f1
        best_loss_thr <- loss_scores[i]
        best_gain_thr <- gain_scores[j]
      }
    }
  }

  pred_class_best <- ifelse(pred_score < best_loss_thr, -1L,
                            ifelse(pred_score > best_gain_thr, 1L, 0L))

  tp_g <- sum(true_class == 1L & pred_class_best == 1L)
  sens_g <- if (sum(true_class == 1L) > 0) tp_g / sum(true_class == 1L) else NA_real_
  prec_g <- if (sum(pred_class_best == 1L) > 0) tp_g / sum(pred_class_best == 1L) else NA_real_

  tp_l <- sum(true_class == -1L & pred_class_best == -1L)
  sens_l <- if (sum(true_class == -1L) > 0) tp_l / sum(true_class == -1L) else NA_real_
  prec_l <- if (sum(pred_class_best == -1L) > 0) tp_l / sum(pred_class_best == -1L) else NA_real_

  c(max_f1 = best_f1,
    gain_thr = best_gain_thr, sens_gain = sens_g, prec_gain = prec_g,
    loss_thr = best_loss_thr, sens_loss = sens_l, prec_loss = prec_l)
}

# Bin-level pooled metrics. CHANGED vs 06: gain/loss truth labels use a CN band
# (gain_cn_thr / loss_cn_thr) instead of the hard >2 / <2 cut.
calc_metrics <- function(pred_mat, true_mat_int, tumor_ids, EPS, max_candidates = 200,
                         gain_cn_thr = 2.25, loss_cn_thr = 1.75,
                         metrics = c("correlation","auc_gain","auc_loss","max_f1")) {

  empty_result <- data.frame(
    pearson = NA_real_, spearman = NA_real_, max_f1 = NA_real_,
    auc_gain = NA_real_, pauc_gain = NA_real_,
    auc_loss = NA_real_, pauc_loss = NA_real_
  )

  common_cells <- intersect(rownames(pred_mat), intersect(rownames(true_mat_int), tumor_ids))
  common_genes <- intersect(colnames(pred_mat), colnames(true_mat_int))

  if (length(common_genes) == 0 || length(common_cells) == 0) return(empty_result)

  P_sub <- pred_mat[common_cells,common_genes, drop=FALSE]
  T_sub <- true_mat_int[common_cells,common_genes, drop=FALSE]
  T_sub_log2 <- log2((pmax(T_sub, 0) + EPS) / 2)

  p_vec <- as.numeric(P_sub)
  t_vec <- as.numeric(T_sub)
  t_log2_vec <- as.numeric(T_sub_log2)

  valid <- is.finite(p_vec) & is.finite(t_vec) & is.finite(t_log2_vec)
  p_vec <- p_vec[valid]; t_vec <- t_vec[valid]; t_log2_vec <- t_log2_vec[valid]

  if (length(p_vec) < 50) return(empty_result)

  result <- empty_result

  if ("correlation" %in% metrics) {
    result$pearson <- if (sd(p_vec) > 0 && sd(t_log2_vec) > 0) cor(p_vec, t_log2_vec) else 0
    result$spearman <- if (length(unique(p_vec)) > 1 && length(unique(t_log2_vec)) > 1) {
      cor(p_vec, t_log2_vec, method = "spearman", use = "complete.obs")
    } else 0
  }

  needs_class <- any(c("auc_gain","auc_loss","max_f1") %in% metrics)
  if (needs_class) {
    gt_gain <- t_vec > gain_cn_thr
    gt_loss <- t_vec < loss_cn_thr
    t_class <- rep(0L, length(t_vec))
    t_class[gt_gain] <- 1L
    t_class[gt_loss] <- -1L
  }

  calc_roc_stats <- function(truth_binary, score_vec, score_threshold, direction_gain = TRUE) {
    if (length(unique(truth_binary)) < 2) return(c(auc=NA_real_, pauc=NA_real_))
    if (sum(truth_binary == 1, na.rm=TRUE) == 0) return(c(auc=NA_real_, pauc=NA_real_))
    scores <- if (direction_gain) score_vec else -score_vec
    thr <- if (direction_gain) score_threshold else -score_threshold
    keep <- !is.na(scores) & !is.na(truth_binary)
    scores <- scores[keep]; y <- truth_binary[keep]
    roc_obj <- suppressMessages(pROC::roc(response = y, predictor = scores, quiet=TRUE))
    auc_val <- as.numeric(pROC::auc(roc_obj))
    pred_pos <- ifelse(scores > thr, 1, 0)
    s_max <- sum(pred_pos==1 & y==1) / sum(y==1)
    pauc_val <- if (!is.na(s_max) && s_max > 0) {
      as.numeric(pROC::auc(roc_obj, partial.auc = c(0, s_max),
                           partial.auc.focus = "sens", partial.auc.correct = FALSE))
    } else 0
    c(auc=auc_val, pauc=pauc_val)
  }

  if ("auc_gain" %in% metrics) {
    stats_gain <- calc_roc_stats(as.numeric(gt_gain), p_vec, score_threshold=0, direction_gain = TRUE)
    result$auc_gain <- stats_gain["auc"]
    result$pauc_gain <- stats_gain["pauc"]
  }
  if ("auc_loss" %in% metrics) {
    stats_loss <- calc_roc_stats(as.numeric(gt_loss), p_vec, score_threshold=0, direction_gain = FALSE)
    result$auc_loss <- stats_loss["auc"]
    result$pauc_loss <- stats_loss["pauc"]
  }
  if ("max_f1" %in% metrics) {
    opt_f1 <- evaluate_optimal_f1_score_log2ratio(
      pred_score = p_vec, true_class = t_class,
      baseline = 0, max_candidates = max_candidates
    )
    result$max_f1 <- opt_f1["max_f1"]
  }

  result
}

# NEW: per-spot Spearman/Pearson of predicted vs true-log2 profile (across genes),
# plus per-spot event counts under the CN band. One row per spot.
calc_perspot_metrics <- function(pred_mat, gt_mat, meta_df, EPS,
                                 gain_cn_thr = 2.25, loss_cn_thr = 1.75,
                                 min_genes = 20) {
  spots <- rownames(pred_mat)
  meta_idx <- match(spots, meta_df$spot_id)

  gt_log2 <- log2((pmax(gt_mat, 0) + EPS) / 2)

  n <- length(spots)
  spearman <- rep(NA_real_, n)
  pearson  <- rep(NA_real_, n)
  n_eval   <- integer(n)
  n_gain   <- integer(n)
  n_loss   <- integer(n)

  for (i in seq_len(n)) {
    p  <- pred_mat[i, ]
    cn <- gt_mat[i, ]
    tl <- gt_log2[i, ]
    ok <- is.finite(p) & is.finite(tl) & is.finite(cn)
    ne <- sum(ok)
    n_eval[i] <- ne
    if (ne < min_genes) next
    po <- p[ok]; tlo <- tl[ok]; cno <- cn[ok]
    n_gain[i] <- sum(cno > gain_cn_thr)
    n_loss[i] <- sum(cno < loss_cn_thr)
    if (sd(po) > 0 && sd(tlo) > 0) pearson[i] <- cor(po, tlo)
    if (length(unique(po)) > 1 && length(unique(tlo)) > 1) {
      spearman[i] <- cor(po, tlo, method = "spearman")
    }
  }

  data.frame(
    spot_id = spots,
    tumor_purity_fraction = meta_df$tumor_purity_fraction[meta_idx],
    purity_bin = meta_df$purity_bin[meta_idx],
    spearman = spearman,
    pearson = pearson,
    n_eval_genes = n_eval,
    n_gain_genes = n_gain,
    n_loss_genes = n_loss,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--root_dir"), type = "character", default = "<DATA_ROOT>/synthetic/processed_data/random_batch_BAF/generate_code/Complete_Batch_v2", help = "Root directory containing Rep_* folders"),
  make_option(c("--out_dir"), type = "character", default = "<DATA_ROOT>/synthetic/results/task1_cnv_profile/purity_gradient_eval_06b", help = "Output directory"),
  make_option(c("--min_purity"), type = "double", default = 0.2, help = "Min purity for bins"),
  make_option(c("--max_purity"), type = "double", default = 1.0, help = "Max purity for bins"),
  make_option(c("--bin_width"), type = "double", default = 0.1, help = "Bin width"),
  make_option(c("--gain_cn_thr"), type = "double", default = 2.25,
              help = "Diluted-CN above this = TRUE gain (detectability band) [default %default]"),
  make_option(c("--loss_cn_thr"), type = "double", default = 1.75,
              help = "Diluted-CN below this = TRUE loss (detectability band) [default %default]"),
  make_option(c("--min_genes_perspot"), type = "integer", default = 20,
              help = "Min finite genes required to compute a per-spot metric [default %default]"),
  make_option(c("--max_candidates"), type = "integer", default = 200,
              help = "Max number of gain/loss threshold candidates for optimal-F1 grid search [default %default]"),
  make_option(c("--tool_out_name"), type = "character", default = "Tool_Outputs",
              help = "Tool output folder name under each rep [default %default]"),
  make_option(c("--std_in_name"), type = "character", default = "Standardized_Inputs",
              help = "Standardized input folder name under each rep [default %default]"),
  make_option(c("--metrics"), type = "character", default = "correlation,auc_gain,auc_loss,max_f1",
              help = "Comma-separated bin-level metrics. Choices: correlation,auc_gain,auc_loss,max_f1 [default %default]"),
  make_option(c("--tools"), type = "character", default = NULL,
              help = "Comma-separated tool folder names to evaluate. Default: all tools found.")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$root_dir)) stop("Usage: Rscript evaluate_purity_gradient.R --root_dir <PATH>")
if (!dir.exists(opt$root_dir)) stop("root_dir does not exist: ", opt$root_dir)
if (!is.finite(opt$max_candidates) || opt$max_candidates < 2) {
  stop("--max_candidates must be an integer >= 2 (got ", opt$max_candidates, ")")
}
if (!(opt$loss_cn_thr < 2 && opt$gain_cn_thr > 2)) {
  stop(sprintf("--loss_cn_thr (%.3f) must be <2 and --gain_cn_thr (%.3f) must be >2", opt$loss_cn_thr, opt$gain_cn_thr))
}
message(sprintf("[Config] CN band: loss < %.2f, neutral, gain > %.2f", opt$loss_cn_thr, opt$gain_cn_thr))
message(sprintf("[Config] Optimal-F1 grid-search max_candidates = %d", opt$max_candidates))

ALL_METRICS <- c("correlation", "auc_gain", "auc_loss", "max_f1")
selected_metrics <- trimws(strsplit(opt$metrics, ",", fixed = TRUE)[[1]])
selected_metrics <- selected_metrics[nzchar(selected_metrics)]
unknown_metrics <- setdiff(selected_metrics, ALL_METRICS)
if (length(unknown_metrics) > 0) stop("Unknown metric(s) in --metrics: ", paste(unknown_metrics, collapse = ","))
if (length(selected_metrics) == 0) stop("--metrics produced empty list")
message(sprintf("[Config] Bin-level metrics: %s", paste(selected_metrics, collapse = ", ")))

requested_tools <- NULL
if (!is.null(opt$tools) && nzchar(opt$tools)) {
  requested_tools <- trimws(strsplit(opt$tools, ",", fixed = TRUE)[[1]])
  requested_tools <- requested_tools[nzchar(requested_tools)]
}
message(sprintf("[Config] Tools filter: %s",
                if (is.null(requested_tools)) "<all>" else paste(requested_tools, collapse = ", ")))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Build per-tool global Spot x Gene matrices (same loop as 06)
# ------------------------------------------------------------------------------
rep_dirs <- list.dirs(opt$root_dir, full.names = TRUE, recursive = FALSE)
rep_dirs <- rep_dirs[grepl("^Rep_", basename(rep_dirs))]
message(sprintf("[Init] Found %d Rep_* folders.", length(rep_dirs)))

global_meta_all <- list()
global_pred_by_tool <- list()
global_gt_by_tool <- list()
global_meta_by_tool <- list()

for (rep_dir in rep_dirs) {
  sample_id <- basename(rep_dir)

  meta_file <- file.path(rep_dir, opt$std_in_name, "Metadata_cellprop.csv")
  gt_file <- file.path(rep_dir, "synthetic", "spatial_data", "st-genome_profile.parquet")
  tool_root <- file.path(rep_dir, opt$tool_out_name)

  if (!file.exists(meta_file) || !file.exists(gt_file) || !dir.exists(tool_root)) {
    message(sprintf("  [%s] Skip (missing metadata/gt/tool outputs).", sample_id))
    next
  }

  meta_df <- fread(meta_file)
  purity_col <- 'Tumor_prop'
  spot_col <- colnames(meta_df)[1]
  meta_rep <- data.frame(
    sample_id = sample_id,
    spot_local = as.character(meta_df[[spot_col]]),
    tumor_purity_fraction = suppressWarnings(as.numeric(meta_df[[purity_col]])),
    stringsAsFactors = FALSE
  )
  meta_rep$spot_id <- paste(meta_rep$sample_id, meta_rep$spot_local, sep = "_")
  meta_rep$purity_bin <- make_purity_bin(
    meta_rep$tumor_purity_fraction,
    min_p = opt$min_purity, max_p = opt$max_purity, by = opt$bin_width
  )
  global_meta_all[[sample_id]] <- meta_rep

  gt_mat <- read_gt_matrix(gt_file)
  if (is.null(rownames(gt_mat))) {
    message(sprintf("  [%s] Skip (GT matrix has no rownames).", sample_id))
    next
  }
  rownames(gt_mat) <- paste(sample_id, rownames(gt_mat), sep = "_")

  tools <- list.dirs(tool_root, full.names = FALSE, recursive = FALSE)
  if (!is.null(requested_tools)) tools <- intersect(tools, requested_tools)

  for (tool in tools) {
    tool_out_dir <- file.path(tool_root, tool)
    eval_specs <- get_tool_eval_specs(tool, tool_out_dir)
    if (length(eval_specs) == 0) next

    for (spec in eval_specs) {
      if (has_global_cache(opt$out_dir, spec$label)) next

      pred_mat <- read_prediction_matrix(spec$file, spec$matrix_type)
      if (is.null(pred_mat)) next
      rownames(pred_mat) <- paste(sample_id, rownames(pred_mat), sep = "_")

      common_spots <- intersect(rownames(pred_mat), intersect(rownames(gt_mat), meta_rep$spot_id))
      common_genes <- intersect(colnames(pred_mat), colnames(gt_mat))
      if (length(common_spots) < 5 || length(common_genes) == 0) next

      pred_sub <- pred_mat[common_spots, common_genes, drop = FALSE]
      gt_sub <- gt_mat[common_spots, common_genes, drop = FALSE]
      meta_sub <- meta_rep[match(common_spots, meta_rep$spot_id), , drop = FALSE]

      global_pred_by_tool[[spec$label]] <- rbind_union_cols(global_pred_by_tool[[spec$label]], pred_sub)
      global_gt_by_tool[[spec$label]] <- rbind_union_cols(global_gt_by_tool[[spec$label]], gt_sub)
      global_meta_by_tool[[spec$label]] <- bind_rows(global_meta_by_tool[[spec$label]], meta_sub)
    }
  }
}

if (length(global_meta_all) > 0) {
  meta_all_df <- bind_rows(global_meta_all)
  fwrite(meta_all_df, file.path(opt$out_dir, "global_metadata_all_spots.csv"))
}

# ------------------------------------------------------------------------------
# Per-tool: write caches, bin-level pooled metrics, AND per-spot metrics
# ------------------------------------------------------------------------------
bin_metrics_all <- list()
perspot_all <- list()
tool_summary <- list()

cached_tool_dirs <- list.dirs(opt$out_dir, recursive = FALSE, full.names = FALSE)
cached_tools <- cached_tool_dirs[vapply(cached_tool_dirs, function(t) has_global_cache(opt$out_dir, t), logical(1))]
cached_tools <- setdiff(cached_tools, SPLIT_TOOL_NAMES)
all_tools <- union(names(global_pred_by_tool), cached_tools)

if (!is.null(requested_tools)) {
  all_tools <- all_tools[vapply(all_tools, function(t) {
    base <- sub("_(CNV|Expr)$", "", t)
    base %in% requested_tools || t %in% requested_tools
  }, logical(1))]
}

for (tool in all_tools) {
  use_cache <- has_global_cache(opt$out_dir, tool)

  if (use_cache) {
    cache_obj <- load_global_cache(opt$out_dir, tool)
    pred_global <- cache_obj$pred
    gt_global <- cache_obj$gt
    meta_global <- cache_obj$meta
    message(sprintf("[Cache] %-30s loaded existing global matrices", tool))
  } else {
    pred_global <- global_pred_by_tool[[tool]]
    gt_global <- global_gt_by_tool[[tool]]
    meta_global <- global_meta_by_tool[[tool]]
  }

  if (is.null(pred_global) || is.null(gt_global) || is.null(meta_global)) next

  common_spots <- intersect(rownames(pred_global), rownames(gt_global))
  common_spots <- intersect(common_spots, meta_global$spot_id)
  common_genes <- intersect(colnames(pred_global), colnames(gt_global))
  if (length(common_spots) < 5 || length(common_genes) == 0) next

  pred_global <- pred_global[common_spots, common_genes, drop = FALSE]
  gt_global <- gt_global[common_spots, common_genes, drop = FALSE]
  meta_global <- meta_global[match(common_spots, meta_global$spot_id), , drop = FALSE]

  tool_dir <- file.path(opt$out_dir, tool)
  dir.create(tool_dir, recursive = TRUE, showWarnings = FALSE)

  if (!use_cache) {
    write_matrix_csv(pred_global, file.path(tool_dir, "global_pred_mat.csv.gz"))
    write_matrix_csv(gt_global, file.path(tool_dir, "global_gt_mat.csv.gz"))
    fwrite(meta_global, file.path(tool_dir, "global_metadata.csv"))
  }

  # ---- (1) Per-spot metrics ----
  ps <- calc_perspot_metrics(pred_global, gt_global, meta_global, EPS,
                             gain_cn_thr = opt$gain_cn_thr,
                             loss_cn_thr = opt$loss_cn_thr,
                             min_genes = opt$min_genes_perspot)
  ps$tool <- tool
  fwrite(ps, file.path(tool_dir, "perspot_metrics.csv"))
  perspot_all[[length(perspot_all) + 1]] <- ps

  # ---- (2) Bin-level pooled metrics (band-aware) ----
  this_tool_metrics <- list()
  bin_labels <- sort(unique(meta_global$purity_bin))
  bin_labels <- bin_labels[!is.na(bin_labels) & bin_labels != ""]
  for (b in bin_labels) {
    ids <- meta_global %>% dplyr::filter(purity_bin == b) %>% dplyr::pull(spot_id)
    bin_ids <- intersect(ids, rownames(pred_global))
    if (length(bin_ids) < 5) next

    pred_bin <- pred_global[bin_ids, , drop = FALSE]
    gt_bin <- gt_global[bin_ids, , drop = FALSE]
    m <- calc_metrics(pred_bin, gt_bin, rownames(pred_bin), EPS,
                      max_candidates = opt$max_candidates,
                      gain_cn_thr = opt$gain_cn_thr, loss_cn_thr = opt$loss_cn_thr,
                      metrics = selected_metrics)
    m$tool <- tool
    m$purity_bin <- b
    m$n_spots <- length(bin_ids)
    m$n_genes <- ncol(pred_global)
    this_tool_metrics[[length(this_tool_metrics) + 1]] <- m
    bin_metrics_all[[length(bin_metrics_all) + 1]] <- m
  }
  fwrite(bind_rows(this_tool_metrics), file.path(tool_dir, "purity_bin_metrics.csv"))

  tool_summary[[length(tool_summary) + 1]] <- data.frame(
    tool = tool, n_spots = nrow(pred_global), n_genes = ncol(pred_global),
    n_spots_with_perspot = sum(!is.na(ps$spearman)),
    stringsAsFactors = FALSE
  )
  message(sprintf("[Done] %-30s spots=%d genes=%d perspot_ok=%d",
                  tool, nrow(pred_global), ncol(pred_global), sum(!is.na(ps$spearman))))
}

tool_summary_df <- bind_rows(tool_summary)
bin_metrics_df <- bind_rows(bin_metrics_all)
perspot_df <- bind_rows(perspot_all)

fwrite(tool_summary_df, file.path(opt$out_dir, "tool_matrix_summary.csv"))
fwrite(bin_metrics_df, file.path(opt$out_dir, "purity_bin_metrics_all_tools.csv"))
fwrite(perspot_df, file.path(opt$out_dir, "perspot_metrics_all_tools.csv"))

# ------------------------------------------------------------------------------
# Plot 1: per-spot Spearman vs CONTINUOUS purity (scatter + loess), faceted by tool
# ------------------------------------------------------------------------------
if (nrow(perspot_df) > 0 && any(!is.na(perspot_df$spearman))) {
  pp <- perspot_df %>% dplyr::filter(is.finite(spearman), is.finite(tumor_purity_fraction))
  p_ps <- ggplot(pp, aes(x = tumor_purity_fraction, y = spearman)) +
    geom_point(alpha = 0.15, size = 0.5, colour = "grey40") +
    geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.9) +
    facet_wrap(~ tool, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(
      title = "Per-spot CNV accuracy vs tumour purity (approach A: vs diluted GT)",
      subtitle = "Each point = one spot's Spearman(pred, true-log2) across genes; red = loess",
      x = "Tumour purity fraction (per spot)",
      y = "Per-spot Spearman"
    ) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(opt$out_dir, "perspot_spearman_vs_purity.png"), p_ps,
         width = 14, height = 9, dpi = 300)
  message("[Done] Saved plot: perspot_spearman_vs_purity.png")
} else {
  message("[Skip] No finite per-spot Spearman values to plot.")
}

# ------------------------------------------------------------------------------
# Plot 2: bin-level pooled metrics vs purity bin (same as 06)
# ------------------------------------------------------------------------------
metric_plot_cols <- intersect(c("spearman", "max_f1"), colnames(bin_metrics_df))
if (nrow(bin_metrics_df) > 0 && length(metric_plot_cols) > 0) {
  parse_bin_start <- function(x) suppressWarnings(as.numeric(sub("-.*$", "", x)))
  bin_levels <- unique(as.character(bin_metrics_df$purity_bin))
  bin_starts <- parse_bin_start(bin_levels)
  bin_levels <- if (all(is.finite(bin_starts))) bin_levels[order(bin_starts)] else sort(bin_levels)

  plot_data <- bin_metrics_df %>%
    mutate(purity_bin = factor(as.character(purity_bin), levels = bin_levels)) %>%
    pivot_longer(cols = all_of(metric_plot_cols), names_to = "metric", values_to = "score")

  p_purity <- ggplot(plot_data, aes(x = purity_bin, y = score, color = tool, group = tool)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(aes(size = n_spots), alpha = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_size_continuous(range = c(1, 5), name = "Num Spots") +
    theme_bw(base_size = 12) +
    labs(
      title = "Tool performance vs tumour purity bin (pooled, band-aware truth)",
      subtitle = sprintf("Truth band: loss<%.2f, gain>%.2f CN; dot size = #spots", opt$loss_cn_thr, opt$gain_cn_thr),
      x = "Tumour purity bin", y = "Performance score"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(opt$out_dir, "purity_bin_metrics_trend.png"), p_purity, width = 12, height = 8, dpi = 300)
  message("[Done] Saved plot: purity_bin_metrics_trend.png")
} else {
  message("[Skip] Not enough data to draw purity trend plot.")
}

message(sprintf("\nFinished. Outputs saved to: %s", opt$out_dir))
