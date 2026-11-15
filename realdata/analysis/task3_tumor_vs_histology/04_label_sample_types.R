#!/usr/bin/env Rscript

# Add Task-3 tumor-composition labels to the combined metrics table.
#
# Sample composition is derived from the least-filtered run, AllSpots_no-ref,
# using max(n_truth_tumor) and max(n_truth_normal) across tools.
# Input and outputs are fixed to:
#   realdata_unified_by_dataset/results/task3_tumor_vs_histology/metrics_combined.csv
#   realdata_unified_by_dataset/results/task3_tumor_vs_histology/metrics_combined_relabelled.csv
#   realdata_unified_by_dataset/results/task3_tumor_vs_histology/sample_type_summary.csv

suppressPackageStartupMessages({
  library(data.table)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

realdata_root <- normalizePath(file.path(script_dir, "..", "..", ".."),
                               winslash = "/", mustWork = TRUE)
result_dir <- file.path(realdata_root, "results", "task3_tumor_vs_histology")
input_csv <- file.path(result_dir, "metrics_combined.csv")
output_csv <- file.path(result_dir, "metrics_combined_relabelled.csv")
sample_summary_csv <- file.path(result_dir, "sample_type_summary.csv")

composition_run <- "AllSpots_no-ref"
balanced_lower <- 0.10
balanced_upper <- 0.90

assign_sample_type <- function(n_tumor,
                               n_normal,
                               balanced_lower = 0.10,
                               balanced_upper = 0.90) {
  n_tumor <- as.numeric(n_tumor)
  n_normal <- as.numeric(n_normal)
  total <- n_tumor + n_normal
  tumor_proportion <- ifelse(total > 0, n_tumor / total, NA_real_)

  out <- rep("other", length(total))
  out[total == 0] <- "no_truth_label"
  out[n_tumor == 0 & n_normal > 0] <- "all_normal"
  out[n_tumor > 0 & n_normal == 0] <- "all_tumor"
  out[tumor_proportion < balanced_lower & n_tumor > 0] <- "normal_dominant"
  out[tumor_proportion >= balanced_lower & tumor_proportion <= balanced_upper] <- "balanced"
  out[tumor_proportion > balanced_upper & n_normal > 0] <- "tumor_dominant"
  out
}

main <- function() {
  if (!file.exists(input_csv)) stop("Input not found: ", input_csv)

  df <- fread(input_csv)
  req <- c("dataset_id", "sample", "run", "tool", "n_truth_tumor", "n_truth_normal")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  message("[Input] ", input_csv)
  message("[Metrics] loaded ",
          uniqueN(paste(df$dataset_id, df$sample, sep = "\r")),
          " samples and ",
          nrow(df), " metric rows")

  numeric_cols <- intersect(c("n_truth_tumor", "n_truth_normal", "n_truth_epi"), names(df))
  for (cc in numeric_cols) df[, (cc) := as.numeric(get(cc))]
  if (!"n_truth_epi" %in% names(df)) {
    df[, n_truth_epi := n_truth_tumor + n_truth_normal]
  }

  composition_runs <- unique(c(composition_run, "AllSpots_no-ref", "Allspots_no-ref"))
  comp <- df[run %in% composition_runs]
  if (!nrow(comp)) {
    stop("No rows found for composition run. Tried: ", paste(composition_runs, collapse = ", "))
  }

  sample_labels <- comp[, .(
    n_truth_tumor_sample = max(n_truth_tumor, na.rm = TRUE),
    n_truth_normal_sample = max(n_truth_normal, na.rm = TRUE),
    n_truth_epi_sample = max(n_truth_epi, na.rm = TRUE)
  ), by = .(dataset_id, sample)]

  sample_labels[, n_truth_total_sample := n_truth_tumor_sample + n_truth_normal_sample]
  sample_labels[, tumor_proportion := ifelse(
    n_truth_total_sample > 0,
    n_truth_tumor_sample / n_truth_total_sample,
    NA_real_
  )]
  sample_labels[, sample_type := assign_sample_type(
    n_tumor = n_truth_tumor_sample,
    n_normal = n_truth_normal_sample,
    balanced_lower = balanced_lower,
    balanced_upper = balanced_upper
  )]

  drop_cols <- intersect(
    names(df),
    c("sample_type", "tumor_proportion", "n_truth_tumor_sample",
      "n_truth_normal_sample", "n_truth_epi_sample", "n_truth_total_sample",
      "false_positive_rate", "false_negative_rate")
  )
  if (length(drop_cols)) df[, (drop_cols) := NULL]

  df_relabelled <- merge(df, sample_labels, by = c("dataset_id", "sample"), all.x = TRUE)
  if ("specificity" %in% names(df_relabelled)) {
    df_relabelled[, false_positive_rate := ifelse(!is.na(specificity), 1 - specificity, NA_real_)]
  }
  if ("sensitivity" %in% names(df_relabelled)) {
    df_relabelled[, false_negative_rate := ifelse(!is.na(sensitivity), 1 - sensitivity, NA_real_)]
  }

  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(df_relabelled, output_csv)
  fwrite(sample_labels, sample_summary_csv)

  message("[Saved] ", output_csv)
  message("[Saved] ", sample_summary_csv)
  message("[Summary] sample types:")
  print(sample_labels[, .N, by = sample_type][order(sample_type)])
}

main()
