#!/usr/bin/env Rscript
# Combine the per-(sample,run,tool) epithelial metrics into the pooled long-form
# table metrics_combined.csv consumed by 04_label_sample_types.R. This is the
# compute-only half of the former 05_plot_performance_combined.R (the side
# artefact it wrote just before finishing); all plotting has been removed.
#
# Default input:
#   results/task3_tumor_vs_histology/metrics_epithelial/
#       metrics_by_sample_run_tool_epithelial.csv
#   (written by 02_evaluate_epithelial_intersection.R)
# Default output:
#   results/task3_tumor_vs_histology/metrics_combined.csv
#
# Root/overview discovery matches 03 (via task3_sample_overview_helpers.R); pass
# --realdata_root / --sample_overview explicitly in a code-only checkout.
#
# Usage:
#   Rscript 03_combine_metrics.R [--input_csvs a.csv,b.csv] [--dataset_roots r1,r2] \
#           [--realdata_root DIR] [--out_dir DIR] [--sample_overview CSV] \
#           [--include_runs a,b] [--exclude_runs a,b]

suppressPackageStartupMessages(library(data.table))

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

pick_col <- function(dt, candidates, label) {
  hit <- candidates[candidates %in% names(dt)]
  if (!length(hit)) stop("Missing metric column for ", label,
                         ". Tried: ", paste(candidates, collapse = ", "))
  hit[1]
}

resolve_metric_csv <- function(dataset_root) {
  file.path(dataset_root, "evaluation", "task3_tumor_vs_histology",
            "metrics_epithelial", "metrics_by_sample_run_tool_epithelial.csv")
}

main <- function() {
  dataset_roots <- NULL
  in_csvs <- NULL
  realdata_root <- NULL
  out_dir <- NULL
  exclude_runs <- NULL
  include_runs <- NULL
  sample_overview <- NULL

  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("--dataset_roots", "--dataset_root") && i + 1L <= length(args)) {
      dataset_roots <- c(dataset_roots, parse_csv_arg(args[i + 1L])); i <- i + 2L; next
    }
    if (key %in% c("--input_csvs", "--input_csv") && i + 1L <= length(args)) {
      in_csvs <- c(in_csvs, parse_csv_arg(args[i + 1L])); i <- i + 2L; next
    }
    if (key == "--realdata_root" && i + 1L <= length(args)) {
      realdata_root <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--out_dir" && i + 1L <= length(args)) {
      out_dir <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--exclude_runs" && i + 1L <= length(args)) {
      exclude_runs <- c(exclude_runs, parse_csv_arg(args[i + 1L])); i <- i + 2L; next
    }
    if (key == "--include_runs" && i + 1L <= length(args)) {
      include_runs <- c(include_runs, parse_csv_arg(args[i + 1L])); i <- i + 2L; next
    }
    if (key == "--sample_overview" && i + 1L <= length(args)) {
      sample_overview <- args[i + 1L]; i <- i + 2L; next
    }
    i <- i + 1L
  }

  sample_overview <- resolve_sample_overview(
    sample_overview = sample_overview,
    dataset_root = if (!is.null(dataset_roots)) dataset_roots[1] else NULL,
    start_dir = script_dir
  )
  realdata_root <- resolve_realdata_root(
    realdata_root = realdata_root,
    sample_overview = sample_overview,
    start_dir = script_dir
  )
  results_dir <- task3_results_dir(realdata_root)

  input_roots <- NULL
  if (is.null(in_csvs) && !is.null(dataset_roots)) {
    in_csvs <- vapply(dataset_roots, resolve_metric_csv, character(1))
    input_roots <- dataset_roots
  }
  if (is.null(in_csvs)) {
    in_csvs <- file.path(results_dir, "metrics_epithelial",
                         "metrics_by_sample_run_tool_epithelial.csv")
  }
  if (is.null(in_csvs) || !length(in_csvs)) {
    stop("Provide --input_csvs or run 02_evaluate_epithelial_intersection.R first.")
  }
  if (is.null(input_roots)) {
    input_roots <- rep(NA_character_, length(in_csvs))
  }
  if (is.null(out_dir)) out_dir <- results_dir

  dt_list <- list()
  for (idx in seq_along(in_csvs)) {
    p <- in_csvs[idx]
    if (!file.exists(p)) {
      message("[Warn] missing metrics CSV: ", p); next
    }
    one <- fread(p)
    if (!"dataset_id" %in% names(one)) {
      one[, dataset_id := basename(dirname(dirname(dirname(dirname(p)))))]
    }
    one[, sample := normalize_sample_id(sample)]

    root <- input_roots[idx]
    if (!is.na(root) && nzchar(root)) {
      overview_dataset_id <- basename(normalizePath(root, winslash = "/", mustWork = FALSE))
      one <- filter_dt_to_histology_samples(
        dt = one, sample_overview = sample_overview, dataset_id = overview_dataset_id)
    } else {
      one <- filter_dt_to_histology_samples(
        dt = one, sample_overview = sample_overview, dataset_id = NULL)
    }
    dt_list[[length(dt_list) + 1L]] <- one
  }
  if (!length(dt_list)) stop("No metrics rows loaded.")

  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  req <- c("dataset_id", "sample", "run", "tool")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  if (!is.null(include_runs) && length(include_runs)) {
    dt <- dt[run %in% include_runs]
  }
  if (!is.null(exclude_runs) && length(exclude_runs)) {
    dt <- dt[!(run %in% exclude_runs)]
  }
  if (!nrow(dt)) stop("No metrics rows remain after run filters.")

  # Metric value columns (same names the historical combined table carried).
  col_f1   <- pick_col(dt, c("f1_score", "f1score"), "f1 score")
  col_sens <- pick_col(dt, c("sensitivity"), "sensitivity")
  col_spec <- pick_col(dt, c("specificity", "specifity"), "specificity")
  col_acc  <- pick_col(dt, c("accuracy"), "accuracy")
  col_cov  <- if ("coverage_epi" %in% names(dt)) "coverage_epi" else NA_character_

  dt[, f1_value := as.numeric(get(col_f1))]
  dt[, sensitivity_value := as.numeric(get(col_sens))]
  dt[, specificity_value := as.numeric(get(col_spec))]
  dt[, accuracy_value := as.numeric(get(col_acc))]
  if (!is.na(col_cov)) dt[, coverage_value := as.numeric(get(col_cov))]
  dt[, sample_global := paste0(dataset_id, ":", sample)]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(out_dir, "metrics_combined.csv")
  fwrite(dt, out_csv)
  message("[Saved] ", out_csv, " (", nrow(dt), " rows)")
}

main()
