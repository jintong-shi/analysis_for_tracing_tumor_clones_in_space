#!/usr/bin/env Rscript

# Evaluate tumor vs histology on epithelial spots only.
#
# The evaluation universe is fixed to histology-defined *epithelial* spots
# (truth_label in {Tumor, Normal}), then intersected with spots the tool produced an
# inferred label for (pred_bin in {0, 1}). Spots flagged Others/stroma/etc. in
# histology, or returned as not.defined/unknown/Unassigned by a tool, are
# dropped from the denominator rather than counted as misses.

suppressPackageStartupMessages({
  library(data.table)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "task3_sample_overview_helpers.R"))

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (prec + rec) == 0) return(NA_real_)
  (1 + beta^2) * (prec * rec) / (beta^2 * prec + rec)
}

# Computes metrics on the epithelial-and-inferred intersection.
# n_truth_epi   : truth has truth_bin in {0,1} (epithelial)
# n_pred_inf    : among epithelial spots, tool returned pred_bin in {0,1}
# n_used        : n_pred_inf (the intersection used for metrics)
calc_metrics <- function(dt) {
  epi <- dt[!is.na(truth_bin)]
  use <- epi[!is.na(pred_bin)]
  n_truth_epi <- nrow(epi)
  n_pred_inf <- nrow(use)

  empty <- data.table(
    n_truth_epi = n_truth_epi, n_pred_inf = n_pred_inf, n_used = n_pred_inf,
    coverage_epi = if (n_truth_epi > 0) n_pred_inf / n_truth_epi else NA_real_,
    TP = NA_integer_, TN = NA_integer_, FP = NA_integer_, FN = NA_integer_,
    sensitivity = NA_real_, specificity = NA_real_, precision = NA_real_,
    f1_score = NA_real_, f2_score = NA_real_,
    accuracy = NA_real_, balanced_accuracy = NA_real_,
    n_truth_tumor = 0L, n_truth_normal = 0L
  )
  if (n_pred_inf == 0L) return(empty)

  t <- as.integer(use$truth_bin)
  p <- as.integer(use$pred_bin)

  TP <- sum(p == 1L & t == 1L)
  TN <- sum(p == 0L & t == 0L)
  FP <- sum(p == 1L & t == 0L)
  FN <- sum(p == 0L & t == 1L)

  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  bacc <- mean(c(sens, spec), na.rm = TRUE)

  data.table(
    n_truth_epi = n_truth_epi,
    n_pred_inf = n_pred_inf,
    n_used = n_pred_inf,
    coverage_epi = n_pred_inf / n_truth_epi,
    TP = TP, TN = TN, FP = FP, FN = FN,
    sensitivity = sens,
    specificity = spec,
    precision = prec,
    f1_score = f_beta(prec, sens, 1),
    f2_score = f_beta(prec, sens, 2),
    accuracy = acc,
    balanced_accuracy = bacc,
    n_truth_tumor = sum(t == 1L),
    n_truth_normal = sum(t == 0L)
  )
}

main <- function() {
  dataset_root <- NULL
  realdata_root <- NULL
  dataset_ids <- NULL
  in_csv <- NULL
  out_dir <- NULL
  sample_overview <- NULL

  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--dataset_root" && i + 1L <= length(args)) {
      dataset_root <- args[i + 1L]
      dataset_ids <- unique(c(dataset_ids, basename(normalizePath(dataset_root, winslash = "/", mustWork = FALSE))))
      i <- i + 2L; next
    }
    if (key == "--realdata_root" && i + 1L <= length(args)) {
      realdata_root <- args[i + 1L]; i <- i + 2L; next
    }
    if (key %in% c("--dataset_id", "--dataset_ids", "--datasets") && i + 1L <= length(args)) {
      dataset_ids <- unique(c(dataset_ids, parse_csv_arg(args[i + 1L]))); i <- i + 2L; next
    }
    if (key == "--input_csv" && i + 1L <= length(args)) {
      in_csv <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--out_dir" && i + 1L <= length(args)) {
      out_dir <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--sample_overview" && i + 1L <= length(args)) {
      sample_overview <- args[i + 1L]; i <- i + 2L; next
    }
    i <- i + 1L
  }

  sample_overview <- resolve_sample_overview(
    sample_overview = sample_overview,
    dataset_root = dataset_root,
    start_dir = script_dir
  )
  realdata_root <- resolve_realdata_root(
    realdata_root = realdata_root,
    sample_overview = sample_overview,
    start_dir = script_dir
  )
  results_dir <- task3_results_dir(realdata_root)
  if (is.null(in_csv)) {
    in_csv <- file.path(results_dir, "predictions_standardized",
                        "predictions_tumor_normal_long.csv")
  }
  if (is.null(out_dir)) out_dir <- file.path(results_dir, "metrics_epithelial")

  if (!file.exists(in_csv)) stop("Input not found: ", in_csv)
  dt <- fread(in_csv)
  req <- c("dataset_id", "sample", "run", "tool", "truth_bin", "pred_bin")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  if (!is.null(dataset_ids) && length(dataset_ids)) {
    dt <- dt[dataset_id %in% dataset_ids]
  }
  dt[, sample := normalize_sample_id(sample)]
  dt <- filter_dt_to_histology_samples(
    dt = dt,
    sample_overview = sample_overview,
    dataset_id = NULL
  )
  if (!nrow(dt)) {
    stop("No prediction rows remain after Task3 histology-available sample filtering.")
  }
  message("[Sample overview] retained ",
          uniqueN(paste(dt$dataset_id, dt$sample, sep = "\r")),
          " Task3 histology-available samples")

  # The 01 script already inner-joins truth and prediction by spot_id, so dt
  # represents the per-spot intersection of the two. Here we further restrict
  # to histology-epithelial spots and to spots the tool actually classified.
  by_keys <- c("dataset_id", "sample", "run", "tool")
  m <- dt[, calc_metrics(.SD), by = by_keys]

  by_tool <- m[, .(
    n_groups = .N,
    n_groups_with_data = sum(n_used > 0, na.rm = TRUE),
    median_coverage_epi = median(coverage_epi, na.rm = TRUE),
    median_f1 = median(f1_score, na.rm = TRUE),
    mean_f1 = mean(f1_score, na.rm = TRUE),
    median_sensitivity = median(sensitivity, na.rm = TRUE),
    median_specificity = median(specificity, na.rm = TRUE),
    median_balanced_accuracy = median(balanced_accuracy, na.rm = TRUE)
  ), by = .(dataset_id, tool)]

  by_sample <- m[, .(
    n_groups = .N,
    median_coverage_epi = median(coverage_epi, na.rm = TRUE),
    median_f1 = median(f1_score, na.rm = TRUE),
    median_balanced_accuracy = median(balanced_accuracy, na.rm = TRUE)
  ), by = .(dataset_id, sample)]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_detail <- file.path(out_dir, "metrics_by_sample_run_tool_epithelial.csv")
  out_tool   <- file.path(out_dir, "metrics_summary_by_tool_epithelial.csv")
  out_sample <- file.path(out_dir, "metrics_summary_by_sample_epithelial.csv")

  fwrite(m, out_detail)
  fwrite(by_tool, out_tool)
  fwrite(by_sample, out_sample)

  message("[Saved] ", out_detail)
  message("[Saved] ", out_tool)
  message("[Saved] ", out_sample)
}

main()
