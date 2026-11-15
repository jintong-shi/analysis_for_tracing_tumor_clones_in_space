#!/usr/bin/env Rscript
# Score every tool's normalized output against the ground truth and write one
# metrics row per (rep, tool): correlation, max F1, AUC / pAUC for gain and loss.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tibble)
  library(data.table)
  library(purrr)
  library(pROC)
  library(mclust)
  library(optparse)
  library(parallel)
})
EPS <- 1e-12
SPLIT_TOOL_NAMES <- c("InferCNV", "Numbat")
DEPRECATED_TOOL_LABELS <- c("XClone_CNV", "XClone_Expr")
# ==============================================================================
# 1. Metrics
# ==============================================================================
# ------------------------------------------------------------
# Macro-F1 for 3-class labels: -1 (loss), 0 (neutral), +1 (gain)
# ------------------------------------------------------------
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
    if (is.na(f1)) f1 <- 0   # mimic paper: NA -> 0 (e.g., predicts no class at all)
    
    f1s[k] <- f1
  }
  
  mean(f1s)
}

get_f1_score_unweighted_log2ratio <- function(pred_score, thr_loss, thr_gain, true_class) {
  # Discretize
  pred_class <- ifelse(pred_score < thr_loss, -1L,
                       ifelse(pred_score > thr_gain,  1L, 0L))
  
  macro_f1_3class(pred_class = pred_class, true_class = true_class)
}

# ------------------------------------------------------------
# Grid-search optimal thresholds for log2ratio scores
# Returns: best macro-F1, best gain_thr, gain recall/precision, best loss_thr, loss recall/precision
# ------------------------------------------------------------
evaluate_optimal_f1_score_log2ratio <- function(pred_score,
                                                true_class,
                                                baseline = 0,
                                                max_candidates = 200) {
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
  
  # ------------------------------------------------------------
  # Sort once, then cumulative counts by true class.
  # Strictly preserve original rule:
  #   pred_score < loss_thr -> loss
  #   pred_score > gain_thr -> gain
  #   otherwise -> neutral
  # ------------------------------------------------------------
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
  
  # For score < loss_thr:
  # findInterval gives number of unique score groups <= threshold.
  # If threshold exactly equals one group, subtract that group to make it strict <.
  le_loss <- findInterval(loss_scores, u)
  eq_loss <- le_loss > 0L
  eq_loss[eq_loss] <- u[le_loss[eq_loss]] == loss_scores[eq_loss]
  lt_loss <- le_loss - as.integer(eq_loss)
  
  loss_counts <- cum_counts[lt_loss + 1L, , drop = FALSE]
  
  # For score > gain_thr:
  # score > threshold = total - score <= threshold
  le_gain <- findInterval(gain_scores, u)
  gain_counts <- sweep(cum_counts[le_gain + 1L, , drop = FALSE], 2, total_cls, FUN = function(x, t) t - x)
  
  best_f1 <- -Inf
  best_loss_thr <- NA_real_
  best_gain_thr <- NA_real_
  
  # Match expand.grid(loss_thr = loss_scores, gain_thr = gain_scores) tie behavior:
  # gain outer, loss inner.
  for (j in seq_along(gain_scores)) {
    g_thr <- gain_scores[j]
    valid_i <- which(loss_scores <= g_thr)
    if (length(valid_i) == 0) next
    
    gc <- gain_counts[j, ]
    
    for (i in valid_i) {
      lc <- loss_counts[i, ]
      
      # predicted loss
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
      
      # predicted gain
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
      
      # predicted neutral = all - predicted loss - predicted gain
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

# Metrics for one tool on one rep. pred_mat / true_mat_int are Genes x Spots;
# tumor_ids selects the columns to score.
calc_metrics <- function(pred_mat, true_mat_int, tumor_ids, EPS, max_candidates = 200,
                         metrics = c("correlation","auc_gain","auc_loss","max_f1")) {
  
  empty_result <- data.frame(
    pearson = NA_real_, spearman = NA_real_, max_f1 = NA_real_,
    auc_gain = NA_real_, pauc_gain = NA_real_,
    auc_loss = NA_real_, pauc_loss = NA_real_
  )
  
  common_cells <- intersect(rownames(pred_mat), intersect(rownames(true_mat_int),tumor_ids))
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
    gt_gain <- t_vec > 2
    gt_loss <- t_vec < 2
    t_class <- rep(0L, length(t_vec))
    t_class[gt_gain] <- 1L
    t_class[gt_loss] <- -1L
  }
  
  calc_roc_stats <- function(truth_binary, score_vec, score_threshold, direction_gain = TRUE) {
    if (length(unique(truth_binary)) < 2) return(c(auc=NA_real_, pauc=NA_real_))
    if (sum(truth_binary == 1, na.rm=TRUE) == 0) return(c(auc=NA_real_, pauc=NA_real_))
    scores <- if (direction_gain) score_vec else -score_vec
    thr    <- if (direction_gain) score_threshold else -score_threshold
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

state_matrix_to_numeric <- function(state_mat) {
  out <- matrix(0, nrow = nrow(state_mat), ncol = ncol(state_mat),
                dimnames = dimnames(state_mat))
  s <- tolower(as.character(as.matrix(state_mat)))
  out[s %in% c("gain", "amp", "bamp", "hamp", "amplified", "amplification")] <- 1
  out[s %in% c("loss", "del", "bdel", "deletion")] <- -1
  out[s %in% c("neutral", "neu", "loh", "cnloh", "cn_loh")] <- 0
  # Numeric-encoded state files (e.g. STARCH writes -1/0/+1): map by sign so a
  # discrete state matrix is read correctly. Text labels above parse to NA here
  # and are left untouched.
  v <- suppressWarnings(as.numeric(s))
  out[!is.na(v) & v > 0] <- 1
  out[!is.na(v) & v < 0] <- -1
  out[!is.na(v) & v == 0] <- 0
  out
}

read_prediction_matrix <- function(csv_file, matrix_type = c("log2", "state")) {
  matrix_type <- match.arg(matrix_type)
  if (matrix_type == "state") {
    # State CSVs are written via write.csv(..., row.names=TRUE), so use read.csv to preserve headers.
    df_in <- read.csv(csv_file, row.names = 1, check.names = FALSE)
    return(state_matrix_to_numeric(as.matrix(df_in)))
  }
  df_in <- fread(csv_file, data.table = FALSE)
  rownames(df_in) <- df_in[, 1]
  df_in <- df_in[, -1, drop = FALSE]
  as.matrix(df_in)
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

  # STARCH is scored on its discrete state matrix (-1/0/+1), like CalicoST/XClone,
  # rather than on its log2 pseudobulk.
  if (tool == "STARCH" || grepl("^STARCH_", tool)) {
    state_file <- file.path(tool_out_dir, "processed_state_spots_x_genes.csv")
    if (!file.exists(state_file)) return(list())
    return(list(list(label = tool, file = state_file, matrix_type = "state")))
  }

  csv_file <- file.path(tool_out_dir, "processed_log2_spots_x_genes.csv")
  if (!file.exists(csv_file)) return(list())
  list(list(label = tool, file = csv_file, matrix_type = "log2"))
}

read_event_log <- function(f, sep = "\t") {
  if (!file.exists(f)) stop("event_log file not found: ", f)
  ev <- tryCatch(
    read.delim(f, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read event_log: ", e$message)
  )
  if (!"length_mb" %in% names(ev)) {
    if (all(c("Start", "End") %in% names(ev))) {
      ev$length_mb <- (as.numeric(ev$End) - as.numeric(ev$Start)) / 1e6
    } else {
      stop("event_log must contain 'length_mb' or ('Start','End').")
    }
  }
  ev
}

parse_event_log_genes_bp_only <- function(ev_df, cutoff_mb = 3, gene_col = "affected_genes", sep = ";") {
  if (!gene_col %in% names(ev_df)) stop("event_log missing column: ", gene_col)
  if (!"length_mb" %in% names(ev_df)) stop("event_log missing column: length_mb")
  ev_df <- ev_df %>% dplyr::filter(!is.na(.data[[gene_col]]), nzchar(.data[[gene_col]]), is.finite(length_mb))
  
  ev_narrow <- ev_df %>% dplyr::filter(length_mb <= cutoff_mb)
  ev_broad  <- ev_df %>% dplyr::filter(length_mb > cutoff_mb)
  
  split_genes <- function(x) unlist(strsplit(x, sep, fixed = TRUE), use.names = FALSE)
  clean <- function(v) unique(trimws(v[nzchar(v)]))
  
  narrow_genes <- clean(unlist(lapply(ev_narrow[[gene_col]], split_genes)))
  broad_genes  <- clean(unlist(lapply(ev_broad[[gene_col]], split_genes)))
  
  list(
    narrow_genes = narrow_genes,
    broad_genes  = broad_genes,
    cnv_genes    = unique(c(narrow_genes, broad_genes))
  )
}

find_event_log_file <- function(rep_dir) {
  candidates <- c(
    file.path(rep_dir, "event_log.tsv"),
    file.path(rep_dir, "synthetic", "event_log.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

subset_matrix_genes <- function(mat, keep_genes) {
  if (is.null(keep_genes) || length(keep_genes) == 0) return(mat[, 0, drop = FALSE])
  g <- intersect(colnames(mat), keep_genes)
  mat[, g, drop = FALSE]
}

calc_metrics_with_scopes <- function(pred_mat, true_mat, tumor_ids, gene_sets, EPS, max_candidates = 200,
                                     metrics = c("correlation","auc_gain","auc_loss","max_f1")) {
  scopes <- list(
    overall = intersect(colnames(pred_mat), colnames(true_mat)),
    narrow_only = gene_sets$narrow_genes,
    broad_only = gene_sets$broad_genes
  )
  
  out <- lapply(names(scopes), function(scope_name) {
    keep <- scopes[[scope_name]]
    p_sub <- subset_matrix_genes(pred_mat, keep)
    t_sub <- subset_matrix_genes(true_mat, keep)
    n_genes <- length(intersect(colnames(p_sub), colnames(t_sub)))
    m <- calc_metrics(p_sub, t_sub, tumor_ids, EPS, max_candidates = max_candidates, metrics = metrics)
    m$scope <- scope_name
    m$n_genes <- n_genes
    m
  })
  
  bind_rows(out) %>% select(scope, n_genes, everything())
}
# ==============================================================================
# 2. Main Logic
# ==============================================================================
option_list <- list(
  make_option(c("--root_dir"), type="character", default="<DATA_ROOT>/synthetic/processed_data/random_batch_BAF/generate_code/Complete_Batch_v2",help="Root containing Rep_* folders"),
  make_option(c("--out_file"), type="character", default="<DATA_ROOT>/synthetic/results/task1_cnv_profile/benchmark_metrics_normalref.csv", help="Result CSV. Default name reflects the normal-reference run (--tool_out_name Tool_Outputs); override for TME-reference runs"),
  make_option(c("--event_cutoff_mb"), type="double", default=3, help="CNV length cutoff (Mb) for narrow/broad split [default %default]"),
  make_option(c("--tool_out_name"), type="character", default="Tool_Outputs",
              help="Tool output folder name under each rep, e.g. Tool_Outputs or Tool_Output_TME [default %default]"),
  make_option(c("--std_in_name"), type="character", default="Standardized_Inputs",
              help="Standardized input folder name under each rep [default %default]"),
  make_option(c("--max_candidates"), type="integer", default=300,
              help="Max number of gain/loss threshold candidates for optimal-F1 grid search [default %default]"),
  make_option(c("--metrics"), type="character", default="correlation,auc_gain,auc_loss,max_f1",
              help="Comma-separated metrics to compute. Choices: correlation,auc_gain,auc_loss,max_f1 [default %default]"),
  make_option(c("--tools"), type="character", default=NULL,
              help="Comma-separated tool folder names to evaluate (e.g. InferCNV,Numbat). Default: all tools found."),
  make_option(c("--workers"), type="integer", default=4,
              help="Number of Rep-level parallel workers [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$workers) || !is.finite(opt$workers) || opt$workers < 1) {
  stop("--workers must be >= 1")
}
message(sprintf("[Config] Rep-level workers = %d", opt$workers))
if (is.null(opt$root_dir)) stop("Usage: Rscript evaluate_results.R --root_dir <PATH>")
if (!is.null(opt$max_candidates) && (!is.finite(opt$max_candidates) || opt$max_candidates < 2)) {
  stop("--max_candidates must be an integer >= 2 (got ", opt$max_candidates, ")")
}
message(sprintf("[Config] Optimal-F1 grid-search max_candidates = %d", opt$max_candidates))

ALL_METRICS <- c("correlation", "auc_gain", "auc_loss", "max_f1")
selected_metrics <- trimws(strsplit(opt$metrics, ",", fixed = TRUE)[[1]])
selected_metrics <- selected_metrics[nzchar(selected_metrics)]
unknown_metrics <- setdiff(selected_metrics, ALL_METRICS)
if (length(unknown_metrics) > 0) stop("Unknown metric(s) in --metrics: ", paste(unknown_metrics, collapse = ","))
if (length(selected_metrics) == 0) stop("--metrics produced empty list")
message(sprintf("[Config] Metrics enabled: %s", paste(selected_metrics, collapse = ", ")))

requested_tools <- NULL
if (!is.null(opt$tools) && nzchar(opt$tools)) {
  requested_tools <- trimws(strsplit(opt$tools, ",", fixed = TRUE)[[1]])
  requested_tools <- requested_tools[nzchar(requested_tools)]
}
message(sprintf("[Config] Tools filter: %s",
                if (is.null(requested_tools)) "<all>" else paste(requested_tools, collapse = ", ")))

write_checkpoint <- function(df, out_file, label = "checkpoint") {
  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_file <- paste0(out_file, ".tmp")
  write.csv(df, tmp_file, row.names = FALSE)
  ok <- file.rename(tmp_file, out_file)
  if (!ok) {
    file.copy(tmp_file, out_file, overwrite = TRUE)
    unlink(tmp_file)
  }
  message(sprintf("[Checkpoint] Saved %d rows after %s -> %s", nrow(df), label, out_file))
}

reps <- list.dirs(opt$root_dir, full.names = FALSE, recursive = FALSE)
reps <- reps[grep("^Rep_", reps)]
all_reps <- reps
message(sprintf("[Init] Found %d batches to process.", length(all_reps)))

results_df <- data.frame()
completed_pairs <- data.frame(rep = character(0), tool = character(0), stringsAsFactors = FALSE)
if (file.exists(opt$out_file)) {
  checkpoint_df <- tryCatch(fread(opt$out_file, data.table = FALSE), error = function(e) NULL)
  if (!is.null(checkpoint_df) && nrow(checkpoint_df) > 0 && all(c("rep","tool") %in% names(checkpoint_df))) {
    legacy_mask <- checkpoint_df$tool %in% DEPRECATED_TOOL_LABELS
    if (any(legacy_mask)) {
      message(sprintf("[Resume] Dropping %d rows with deprecated tool labels: %s",
                      sum(legacy_mask), paste(DEPRECATED_TOOL_LABELS, collapse = ", ")))
      checkpoint_df <- checkpoint_df[!legacy_mask, , drop = FALSE]
    }
    results_df <- checkpoint_df
    completed_pairs <- unique(checkpoint_df[, c("rep", "tool")])
    completed_pairs$rep <- as.character(completed_pairs$rep)
    completed_pairs$tool <- as.character(completed_pairs$tool)
    message(sprintf("[Resume] Loaded %d rows from %s", nrow(results_df), opt$out_file))
    message(sprintf("[Resume] Completed (rep, tool) pairs in checkpoint: %d", nrow(completed_pairs)))
  } else {
    message(sprintf("[Resume] Existing out_file could not be reused cleanly. Recomputing: %s", opt$out_file))
  }
}

reps <- all_reps

total_reps <- length(reps)
if (total_reps == 0) {
  message("[Resume] No remaining reps to process.")
}

process_one_rep <- function(rep_name,
                            rep_idx,
                            total_reps,
                            opt,
                            requested_tools,
                            completed_pairs,
                            selected_metrics,
                            EPS) {
  rep_dir <- file.path(opt$root_dir, rep_name)
  
  message(sprintf("\n[Rep %d/%d] %s", rep_idx, total_reps, rep_name))
  message(sprintf("  [%s] Loading metadata ...", rep_name))
  
  prop_file <- file.path(rep_dir, opt$std_in_name, "Metadata_cellprop.csv")
  if (!file.exists(prop_file)) {
    alt_files <- c(
      file.path(rep_dir, "Standardized_Inputs", "Metadata_cellprop.csv"),
      file.path(rep_dir, "Standardized_Inputs_TME", "Metadata_cellprop.csv")
    )
    hit <- alt_files[file.exists(alt_files)]
    if (length(hit) > 0) prop_file <- hit[1]
  }
  
  if (!file.exists(prop_file)) {
    message(sprintf("  [%s] Skip: Missing Metadata_cellprop.csv", rep_name))
    return(data.frame())
  }
  
  prop_df <- fread(prop_file)
  
  if (!"Tumor_prop" %in% names(prop_df)) {
    message(sprintf("  [%s] Skip: 'Tumor_prop' column missing in metadata", rep_name))
    return(data.frame())
  }
  
  id_col <- colnames(prop_df)[1]
  target_tumor_spots <- prop_df[Tumor_prop > 0.5][[id_col]]
  
  message(sprintf("  [%s] Target tumor spots (>50%%): %d",
                  rep_name, length(target_tumor_spots)))
  
  gt_file <- file.path(rep_dir, "synthetic", "spatial_data", "st-genome_profile.parquet")
  
  if (!file.exists(gt_file)) {
    message(sprintf("  [%s] Skip: Missing ground truth parquet.", rep_name))
    return(data.frame())
  }
  
  message(sprintf("  [%s] Loading ground truth parquet ...", rep_name))
  gt_df <- read_parquet(gt_file)
  
  if ("index" %in% names(gt_df)) {
    gt_df <- tibble::column_to_rownames(gt_df, "index")
  }
  
  gt_mat <- t(as.matrix(gt_df))
  
  message(sprintf("  [%s] Ground truth matrix ready: %d spots x %d genes",
                  rep_name, nrow(gt_mat), ncol(gt_mat)))
  
  event_file <- find_event_log_file(rep_dir)
  gene_sets <- list(
    narrow_genes = character(0),
    broad_genes = character(0),
    cnv_genes = character(0)
  )
  
  if (!is.null(event_file)) {
    message(sprintf("  [%s] Parsing event log: %s", rep_name, basename(event_file)))
    
    ev_df <- tryCatch(
      read_event_log(event_file, sep = "\t"),
      error = function(e) {
        message(sprintf("  [%s] Warning: failed to parse event_log (%s)",
                        rep_name, e$message))
        NULL
      }
    )
    
    if (!is.null(ev_df)) {
      gene_sets <- parse_event_log_genes_bp_only(
        ev_df,
        cutoff_mb = opt$event_cutoff_mb,
        gene_col = "affected_genes",
        sep = ";"
      )
      
      gene_sets$narrow_genes <- intersect(gene_sets$narrow_genes, colnames(gt_mat))
      gene_sets$broad_genes  <- intersect(gene_sets$broad_genes, colnames(gt_mat))
      gene_sets$cnv_genes    <- intersect(gene_sets$cnv_genes, colnames(gt_mat))
    }
  } else {
    message(sprintf("  [%s] Warning: event_log.tsv not found, narrow/broad scopes will be empty.",
                    rep_name))
  }
  
  message(sprintf("  [%s] Gene scopes ready | CNV: %d | narrow: %d | broad: %d",
                  rep_name,
                  length(gene_sets$cnv_genes),
                  length(gene_sets$narrow_genes),
                  length(gene_sets$broad_genes)))
  
  tool_base_dir <- file.path(rep_dir, opt$tool_out_name)
  
  if (!dir.exists(tool_base_dir)) {
    message(sprintf("  [%s] Skip: Missing %s", rep_name, opt$tool_out_name))
    return(data.frame())
  }
  
  tool_names <- list.dirs(tool_base_dir, full.names = FALSE, recursive = FALSE)
  
  if (!is.null(requested_tools)) {
    tool_names <- intersect(tool_names, requested_tools)
  }
  
  message(sprintf("  [%s] Scanning %d tool folders under %s",
                  rep_name, length(tool_names), opt$tool_out_name))
  
  rep_tool_hits <- 0L
  rep_results_list <- list()
  
  for (tool_idx in seq_along(tool_names)) {
    tool <- tool_names[[tool_idx]]
    tool_out_dir <- file.path(tool_base_dir, tool)
    
    message(sprintf("    [%s | %d/%d] Checking %s",
                    rep_name, tool_idx, length(tool_names), tool))
    
    eval_specs <- get_tool_eval_specs(tool, tool_out_dir)
    
    if (length(eval_specs) == 0) {
      if (tool %in% SPLIT_TOOL_NAMES) {
        message(sprintf(
          "    [%s | %s] Skip: Missing both processed_state_spots_x_genes.csv and processed_log2_spots_x_genes.csv",
          rep_name, tool
        ))
      } else {
        message(sprintf("    [%s | %s] Skip: Missing processed_log2_spots_x_genes.csv",
                        rep_name, tool))
      }
      next
    }
    
    for (spec in eval_specs) {
      if (nrow(completed_pairs) > 0 &&
          any(completed_pairs$rep == rep_name & completed_pairs$tool == spec$label)) {
        message(sprintf("    [%s | %s] Skip: already in checkpoint",
                        rep_name, spec$label))
        next
      }
      
      message(sprintf("    [%s | %s] Reading %s",
                      rep_name, spec$label, basename(spec$file)))
      
      pred_mat <- read_prediction_matrix(spec$file, spec$matrix_type)
      
      message(sprintf("    [%s | %s] Prediction matrix ready: %d spots x %d genes",
                      rep_name, spec$label, nrow(pred_mat), ncol(pred_mat)))
      
      message(sprintf("    [%s | %s] Calculating metrics ...",
                      rep_name, spec$label))
      
      m_df <- calc_metrics_with_scopes(
        pred_mat,
        gt_mat,
        target_tumor_spots,
        gene_sets,
        EPS,
        max_candidates = opt$max_candidates,
        metrics = selected_metrics
      )
      
      m_df$rep <- rep_name
      m_df$tool <- spec$label
      
      rep_results_list[[length(rep_results_list) + 1]] <- m_df
      rep_tool_hits <- rep_tool_hits + 1L
      
      m_overall <- m_df %>% filter(scope == "overall")
      
      if (nrow(m_overall) == 1) {
        message(sprintf("  [%s] Tool: %-14s | Scope: overall | Pearson: %.3f | MaxF1: %.3f",
                        rep_name,
                        spec$label,
                        m_overall$pearson,
                        m_overall$max_f1))
      }
      
      rm(pred_mat)
      gc(verbose = FALSE)
    }
  }
  
  if (length(rep_results_list) > 0) {
    rep_results_df <- bind_rows(rep_results_list)
  } else {
    rep_results_df <- data.frame()
  }
  
  rm(gt_mat, gt_df, prop_df, rep_results_list)
  gc(verbose = FALSE)
  
  message(sprintf("  [%s] Finished: evaluated %d tools with valid processed matrices.",
                  rep_name, rep_tool_hits))
  
  rep_results_df
}

total_reps <- length(reps)

if (total_reps == 0) {
  message("[Resume] No remaining reps to process.")
}

rep_table <- data.frame(
  rep = reps,
  rep_idx = seq_along(reps),
  stringsAsFactors = FALSE
)

workers <- min(opt$workers, max(1L, nrow(rep_table)))

if (workers <= 1L) {
  message("[Parallel] Running sequentially.")
  
  for (ii in seq_len(nrow(rep_table))) {
    one_df <- process_one_rep(
      rep_name = rep_table$rep[ii],
      rep_idx = rep_table$rep_idx[ii],
      total_reps = total_reps,
      opt = opt,
      requested_tools = requested_tools,
      completed_pairs = completed_pairs,
      selected_metrics = selected_metrics,
      EPS = EPS
    )
    
    if (nrow(one_df) > 0) {
      results_df <- bind_rows(results_df, one_df)
      completed_pairs <- unique(results_df[, c("rep", "tool")])
      completed_pairs$rep <- as.character(completed_pairs$rep)
      completed_pairs$tool <- as.character(completed_pairs$tool)
      write_checkpoint(results_df, opt$out_file, label = rep_table$rep[ii])
    } else {
      message(sprintf("  [%s] No valid tool metrics generated; checkpoint unchanged.",
                      rep_table$rep[ii]))
    }
    
    gc(verbose = FALSE)
  }
  
} else {
  message(sprintf("[Parallel] Running Rep-level parallel with %d workers.", workers))
  message("[Parallel] Memory note: only Rep-level is parallelized; tools inside each Rep remain sequential.")
  
  # Chunking prevents launching too many Reps at once and keeps memory controlled.
  chunks <- split(seq_len(nrow(rep_table)), ceiling(seq_len(nrow(rep_table)) / workers))
  
  for (chunk_id in seq_along(chunks)) {
    idxs <- chunks[[chunk_id]]
    
    message(sprintf("[Parallel] Processing chunk %d/%d with %d Rep(s).",
                    chunk_id, length(chunks), length(idxs)))
    
    chunk_results <- parallel::mclapply(
      idxs,
      function(ii) {
        process_one_rep(
          rep_name = rep_table$rep[ii],
          rep_idx = rep_table$rep_idx[ii],
          total_reps = total_reps,
          opt = opt,
          requested_tools = requested_tools,
          completed_pairs = completed_pairs,
          selected_metrics = selected_metrics,
          EPS = EPS
        )
      },
      mc.cores = workers,
      mc.preschedule = FALSE
    )
    
    chunk_df <- bind_rows(chunk_results)
    
    if (nrow(chunk_df) > 0) {
      results_df <- bind_rows(results_df, chunk_df)
      
      completed_pairs <- unique(results_df[, c("rep", "tool")])
      completed_pairs$rep <- as.character(completed_pairs$rep)
      completed_pairs$tool <- as.character(completed_pairs$tool)
      
      write_checkpoint(results_df, opt$out_file, label = paste0("parallel_chunk_", chunk_id))
    } else {
      message(sprintf("[Parallel] Chunk %d generated no valid metrics; checkpoint unchanged.",
                      chunk_id))
    }
    
    rm(chunk_results, chunk_df)
    gc(verbose = FALSE)
  }
}

# ============================================================
# Save
# ============================================================
final_df <- results_df
if (nrow(final_df) > 0) {
  message(sprintf("[Save] Writing %d rows to %s", nrow(final_df), opt$out_file))
  write_checkpoint(final_df, opt$out_file, label = "final save")
  message(sprintf("\n[Done] Results saved to: %s", opt$out_file))
} else {
  message("\n[Done] No valid results generated.")
}
