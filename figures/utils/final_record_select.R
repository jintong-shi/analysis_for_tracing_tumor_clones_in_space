#!/usr/bin/env Rscript
# =====================================================================
# Shared best-case selection for the REAL-data figures (Task 1/2/3).
#
# Centralises the fixed reference-strategy selection (docs/task1_cnv_profile.md)
# that the figures apply to the FULL evaluation tables just before plotting. It
# replaces the former per-task pre-compute scripts: the figure loaders call these
# functions on the full table, so no intermediate *_BEST / *_bestcase CSV exists.
#
# Fixed rule, per (sample, tool):
#   ref-free tools           -> AllSpots_no-ref
#   ref tools with histology -> HistoTumor_Benign-ref  (else Decon_StromaImm-ref)
#   task3 variant            -> every tool -> AllSpots_no-ref, except numbat -> Decon_StromaImm-ref
#
# ASSUMES figure_setup.R was sourced first (final_record_select_task1 uses figure_settings(),
# TOOL_DISPLAY, TOOL_ALIAS_TABLE, read_csv, and the global .figure_setup_dir).
# Source this AFTER figure_setup.R. Defines functions only; no side effects.
#
# sample_overview.csv (sample sheet: sample, datasets, sample_type, bam_available,
# histology_available, ...) is vendored next to this file and is the default for
# the overview-dependent selectors (task1, task2).
# =====================================================================

# Default vendored sample sheet, resolved at call time next to figure_setup.R.
.final_record_select_overview <- function() file.path(.figure_setup_dir, "sample_overview.csv")

# ---- Task 1 ----------------------------------------------------------
# Collapse the full pseudobulk_vs_wes_metrics_all_datasets.csv to one best-case
# row per (dataset, sample, tool). Output shape == historical
# pseudobulk_vs_wes_metrics_BEST.csv (adds dataset_lab/tool_id/tool_lab/sample_type).
final_record_select_task1 <- function(raw, overview = .final_record_select_overview()) {
  METRICS4     <- c("Spearman", "max_F1", "AUC_gain_trunc", "AUC_loss_trunc")
  RUN_ORDER_T1 <- c("AllSpots_no-ref", "Decon_no-ref",
                    "Decon_StromaImm-ref", "HistoTumor_Benign-ref")
  CNV_NEGATIVE <- c("skin_P4_V1", "skin_P4_V2")

  metric_cols <- intersect(METRICS4, colnames(raw))
  if (!setequal(metric_cols, METRICS4))
    stop("final_record_select_task1: missing metric columns: ",
         paste(setdiff(METRICS4, metric_cols), collapse = ", "))

  ds         <- figure_settings("dataset", "real_task1")
  DSID2LABEL <- setNames(ds$label, ds$id)

  alias_lc <- tibble::tibble(key = tolower(TOOL_ALIAS_TABLE$alias),
                             id  = TOOL_ALIAS_TABLE$id) %>%
    distinct(key, .keep_all = TRUE)
  norm_tool_id <- function(x) {
    id <- alias_lc$id[match(tolower(as.character(x)), alias_lc$key)]
    ifelse(is.na(id), as.character(x), id)
  }

  if (!file.exists(overview)) stop("final_record_select_task1: sample_overview not found: ", overview)
  ov <- read_csv(overview, show_col_types = FALSE)
  ov_sample <- trimws(as.character(ov$sample))
  hist_avail <- setNames(
    tolower(trimws(as.character(ov$histology_available))) == "yes", ov_sample)
  sample_type <- setNames(
    ifelse(tolower(trimws(as.character(ov$sample_type))) == "ffpe", "ffpe", "ff"), ov_sample)

  dat0 <- raw %>%
    filter(.data$status == "ok") %>%
    mutate(
      sample      = trimws(as.character(.data$sample)),
      tool_id     = norm_tool_id(.data$tool),
      tool_lab    = unname(TOOL_DISPLAY[.data$tool_id]),
      tool_lab    = ifelse(is.na(.data$tool_lab), as.character(.data$tool), .data$tool_lab),
      dataset_lab = unname(DSID2LABEL[as.character(.data$dataset)]),
      run         = factor(.data$run, levels = RUN_ORDER_T1)
    ) %>%
    filter(!is.na(.data$dataset_lab)) %>%
    filter(.data$sample %in% ov_sample, !.data$sample %in% CNV_NEGATIVE) %>%
    group_by(.data$dataset, .data$dataset_lab, .data$sample, .data$run,
             .data$tool, .data$tool_id, .data$tool_lab) %>%
    summarise(across(all_of(metric_cols),
                     ~ suppressWarnings(max(as.numeric(.x), na.rm = TRUE))),
              .groups = "drop") %>%
    mutate(across(all_of(metric_cols), ~ ifelse(is.infinite(.x), NA_real_, .x)))

  dat0 %>%
    mutate(
      sample_type = unname(sample_type[.data$sample]),
      histo_avail = unname(hist_avail[.data$sample]),
      histo_avail = ifelse(is.na(.data$histo_avail), FALSE, .data$histo_avail),
      run_chr     = as.character(.data$run),
      sel_pri = dplyr::case_when(
        .data$histo_avail  & .data$run_chr == "HistoTumor_Benign-ref" ~ 1L,
        .data$histo_avail  & .data$run_chr == "Decon_StromaImm-ref"   ~ 2L,
        !.data$histo_avail & .data$run_chr == "Decon_StromaImm-ref"   ~ 1L,
        .data$run_chr == "AllSpots_no-ref"                            ~ 3L,
        .data$run_chr == "Decon_no-ref"                               ~ 4L,
        TRUE ~ 9L)
    ) %>%
    group_by(.data$dataset_lab, .data$sample, .data$tool_id, .data$tool_lab) %>%
    slice_min(order_by = .data$sel_pri, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-histo_avail, -run_chr, -sel_pri) %>%
    mutate(run = as.character(.data$run)) %>%
    select(dataset, dataset_lab, sample, sample_type,
           tool, tool_id, tool_lab, run, all_of(METRICS4)) %>%
    arrange(.data$dataset_lab, .data$sample, .data$tool_lab)
}

# ---- Task 2 ----------------------------------------------------------
# Fixed best-case run per (dataset, sample, tool) from the full
# subclone_spatial_assortativity_summary.csv, plus platform group. Output shape
# == historical subclone_spatial_assortativity_bestcase.csv (adds `group`; keeps
# only FFPE + FF-BAM). Runs in the task-2 figure context where data.table is loaded.
final_record_select_task2 <- function(summary, overview = .final_record_select_overview()) {
  REF_FREE_TOOLS <- c("STARCH", "CalicoST", "TUSCAN")
  NEGATIVE_SAMPLES <- data.table::data.table(
    dataset = "skin_tumor_ff", sample = c("skin_P4_V1", "skin_P4_V2"))

  dt <- data.table::as.data.table(summary)
  dt <- dt[!NEGATIVE_SAMPLES, on = c("dataset", "sample")]

  histo <- dt[run == "HistoTumor_Benign-ref",
              .(has_histo = TRUE), by = c("dataset", "sample")]
  dt <- merge(dt, histo, by = c("dataset", "sample"), all.x = TRUE)
  dt[is.na(has_histo), has_histo := FALSE]
  dt[, req_run := data.table::fifelse(
    tool %in% REF_FREE_TOOLS, "AllSpots_no-ref",
    data.table::fifelse(has_histo, "HistoTumor_Benign-ref", "Decon_StromaImm-ref"))]
  best <- dt[run == req_run]
  best[, c("has_histo", "req_run") := NULL]

  if (!file.exists(overview)) stop("final_record_select_task2: sample_overview not found: ", overview)
  ov <- data.table::fread(overview, encoding = "UTF-8")
  data.table::setnames(ov, 1, "sample")
  ov[, sample_type   := tolower(trimws(sample_type))]
  ov[, bam_available := tolower(trimws(bam_available))]
  ov[, group := data.table::fifelse(sample_type == "ffpe", "FFPE",
               data.table::fifelse(bam_available %in% c("yes", "true", "1", "y", "t"),
                                   "FF-BAM", "FF-no-BAM"))]
  grp <- unique(ov[, .(dataset = datasets, sample, group)])
  best <- merge(best, grp, by = c("dataset", "sample"), all.x = TRUE)
  miss <- unique(best[is.na(group), paste(dataset, sample)])
  if (length(miss)) stop("final_record_select_task2: samples missing from sample_overview.csv: ",
                         paste(miss, collapse = ", "))
  best[group %in% c("FFPE", "FF-BAM")]
}

# ---- Task 3 ----------------------------------------------------------
# Best-run-per-tool row filter on metrics_combined_relabelled.csv. Adds no
# columns (identity to the historical *_bestcase.csv).
final_record_select_task3 <- function(d) {
  d %>% filter((.data$tool != "numbat" & .data$run == "AllSpots_no-ref") |
               (.data$tool == "numbat" & .data$run == "Decon_StromaImm-ref"))
}
