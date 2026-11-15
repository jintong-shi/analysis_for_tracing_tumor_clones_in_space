#!/usr/bin/env Rscript

# Standardize CNVkit `.cns` WES segment files (skin_tumor dataset) into the
# unified per-run_sample format consumed by 04_evaluate_pseudobulk_vs_wes_metrics.R.
#
# Input: one CNVkit `.cns` per run_sample. Columns: chromosome, start, end,
# gene, log2, depth, probes, weight (and optionally ci_lo, ci_hi).
#
# Output: <out_dir>/<run_sample>_wes_segments_results.csv with columns
#   run_sample, patient_id, sample, wes_sample_id,
#   CONTIG, START, END, NUM_POINTS_COPY_RATIO, MEAN_LOG2_COPY_RATIO, CALL,
#   depth, weight, gene
#
# CALL is derived from log2 with thresholds: log2 >= gain_threshold -> "+",
# log2 <= loss_threshold -> "-", else "0".

suppressPackageStartupMessages({
  library(data.table)
})

DEFAULT_OUT_BASENAME <- "wes_segments_results.csv"
DEFAULT_GAIN <- 0.3
DEFAULT_LOSS <- -0.3

# Unified-dataset layout under --dataset_dir:
#   <dataset_dir>/samples/sample_table.csv             (run_sample -> cns_path map)
#   <dataset_dir>/truth_links/wes_cnv/                 (default output location)
DATASET_SAMPLE_TABLE <- "samples/sample_table.csv"
DATASET_OUT_SUBDIR   <- "truth_links/wes_cnv"

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

extract_patient_id <- function(run_sample) {
  sub("[-_].*$", "", as.character(run_sample))
}

derive_call <- function(log2, gain_thr, loss_thr) {
  out <- rep("0", length(log2))
  out[is.na(log2)] <- NA_character_
  out[!is.na(log2) & log2 >= gain_thr] <- "+"
  out[!is.na(log2) & log2 <= loss_thr] <- "-"
  out
}

choose_col <- function(nms, candidates, required = TRUE, label = "") {
  nms_key <- normalize_key(nms)
  cnd_key <- normalize_key(candidates)
  idx <- match(cnd_key, nms_key)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    if (required) stop("Cannot find ", label, " column. Tried: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  nms[idx[1]]
}

load_sample_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("sample_table not found: ", path)
  }
  st <- fread(path)
  if (nrow(st) == 0L) stop("sample_table is empty: ", path)

  col_run    <- choose_col(names(st), c("run_sample", "sample", "id"), TRUE, "run_sample")
  col_cns    <- choose_col(names(st), c("cns_path", "cns", "wes_cns", "input_path", "path"), TRUE, "cns_path")
  col_pat    <- choose_col(names(st), c("patient_id", "patient", "case", "case_id"), FALSE)
  col_wesid  <- choose_col(names(st), c("wes_sample_id", "wes_sample", "wes_id", "sample_id"), FALSE)

  st <- st[!is.na(get(col_run)) & nzchar(as.character(get(col_run)))]
  out <- data.table(
    run_sample    = as.character(st[[col_run]]),
    cns_path      = as.character(st[[col_cns]]),
    patient_id    = if (is.na(col_pat))   NA_character_ else as.character(st[[col_pat]]),
    wes_sample_id = if (is.na(col_wesid)) NA_character_ else as.character(st[[col_wesid]])
  )
  out <- out[!is.na(cns_path) & nzchar(cns_path)]
  out
}

read_cns <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE)
  col_chr    <- choose_col(names(dt), c("chromosome", "chr", "chrom", "CONTIG"), TRUE, "chromosome")
  col_start  <- choose_col(names(dt), c("start", "START"), TRUE, "start")
  col_end    <- choose_col(names(dt), c("end", "END"), TRUE, "end")
  col_log2   <- choose_col(names(dt), c("log2", "log2ratio", "MEAN_LOG2_COPY_RATIO"), TRUE, "log2")
  col_probes <- choose_col(names(dt), c("probes", "n_probes", "NUM_POINTS_COPY_RATIO"), FALSE)
  col_depth  <- choose_col(names(dt), c("depth", "mean_depth"), FALSE)
  col_weight <- choose_col(names(dt), c("weight"), FALSE)
  col_gene   <- choose_col(names(dt), c("gene", "genes"), FALSE)

  out <- data.table(
    CONTIG  = as.character(dt[[col_chr]]),
    START   = as.integer(round(as.numeric(dt[[col_start]]))),
    END     = as.integer(round(as.numeric(dt[[col_end]]))),
    MEAN_LOG2_COPY_RATIO = as.numeric(dt[[col_log2]]),
    NUM_POINTS_COPY_RATIO = if (is.na(col_probes)) NA_integer_ else as.integer(round(as.numeric(dt[[col_probes]]))),
    depth   = if (is.na(col_depth))  NA_real_ else as.numeric(dt[[col_depth]]),
    weight  = if (is.na(col_weight)) NA_real_ else as.numeric(dt[[col_weight]]),
    gene    = if (is.na(col_gene))   NA_character_ else as.character(dt[[col_gene]])
  )
  out <- out[!is.na(CONTIG) & !is.na(START) & !is.na(END) & END >= START]
  out
}

log_msg <- function(...) cat(sprintf(...), "\n")

main <- function() {
  sample_table <- NULL
  dataset_dir  <- NULL
  out_dir      <- NULL
  out_basename <- DEFAULT_OUT_BASENAME
  gain_thr     <- DEFAULT_GAIN
  loss_thr     <- DEFAULT_LOSS
  samples      <- NULL
  overwrite    <- FALSE
  dry_run      <- FALSE
  summary_csv  <- NULL

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L) {
    i <- 1L
    while (i <= length(args)) {
      key <- args[i]
      if (key == "--sample_table" && i + 1L <= length(args)) {
        sample_table <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--dataset_dir" && i + 1L <= length(args)) {
        dataset_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_dir" && i + 1L <= length(args)) {
        out_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_basename" && i + 1L <= length(args)) {
        out_basename <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--gain_threshold" && i + 1L <= length(args)) {
        gain_thr <- as.numeric(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--loss_threshold" && i + 1L <= length(args)) {
        loss_thr <- as.numeric(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--samples" && i + 1L <= length(args)) {
        samples <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--summary_csv" && i + 1L <= length(args)) {
        summary_csv <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--overwrite") { overwrite <- TRUE; i <- i + 1L; next }
      if (key == "--dry_run")   { dry_run   <- TRUE; i <- i + 1L; next }
      i <- i + 1L
    }
  }

  if (!is.null(dataset_dir)) {
    if (!dir.exists(dataset_dir)) stop("dataset_dir does not exist: ", dataset_dir)
    if (is.null(sample_table)) sample_table <- file.path(dataset_dir, DATASET_SAMPLE_TABLE)
    if (is.null(out_dir))      out_dir      <- file.path(dataset_dir, DATASET_OUT_SUBDIR)
  }

  if (is.null(sample_table)) stop("--sample_table or --dataset_dir is required")
  if (is.null(out_dir))      stop("--out_dir or --dataset_dir is required")
  if (is.na(gain_thr) || is.na(loss_thr)) stop("gain_threshold/loss_threshold must be numeric")
  if (loss_thr >= gain_thr) stop("loss_threshold must be < gain_threshold")

  st <- load_sample_table(sample_table)
  if (!is.null(samples)) {
    st <- st[normalize_key(run_sample) %in% normalize_key(samples)]
  }
  if (nrow(st) == 0L) stop("no run_samples to process after filtering")

  if (!dir.exists(out_dir) && !dry_run) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  log_msg("[Init] sample_table = %s (%d run_samples)", sample_table, nrow(st))
  log_msg("[Init] out_dir = %s", out_dir)
  log_msg("[Init] thresholds: gain >= %.3f | loss <= %.3f", gain_thr, loss_thr)
  log_msg("[Init] overwrite = %s | dry_run = %s", as.character(overwrite), as.character(dry_run))

  summary <- data.table(
    run_sample = character(), patient_id = character(), wes_sample_id = character(),
    cns_path = character(), n_segments = integer(), n_gain = integer(),
    n_loss = integer(), n_neutral = integer(),
    out_file = character(), status = character()
  )

  for (k in seq_len(nrow(st))) {
    run_sample <- st$run_sample[k]
    cns_path   <- st$cns_path[k]
    patient_id <- if (!is.na(st$patient_id[k]) && nzchar(st$patient_id[k]))
      st$patient_id[k] else extract_patient_id(run_sample)
    wes_id     <- if (!is.na(st$wes_sample_id[k]) && nzchar(st$wes_sample_id[k]))
      st$wes_sample_id[k] else run_sample

    out_file <- file.path(out_dir, paste0(run_sample, "_", out_basename))

    if (!file.exists(cns_path)) {
      log_msg("[Skip] %s | cns_path missing: %s", run_sample, cns_path)
      summary <- rbind(summary, data.table(
        run_sample = run_sample, patient_id = patient_id, wes_sample_id = wes_id,
        cns_path = cns_path, n_segments = 0L, n_gain = 0L, n_loss = 0L, n_neutral = 0L,
        out_file = out_file, status = "missing_cns"
      ), fill = TRUE)
      next
    }

    if (file.exists(out_file) && !overwrite) {
      log_msg("[Skip] %s | output exists: %s", run_sample, out_file)
      summary <- rbind(summary, data.table(
        run_sample = run_sample, patient_id = patient_id, wes_sample_id = wes_id,
        cns_path = cns_path, n_segments = NA_integer_, n_gain = NA_integer_,
        n_loss = NA_integer_, n_neutral = NA_integer_,
        out_file = out_file, status = "exists_skip"
      ), fill = TRUE)
      next
    }

    seg <- read_cns(cns_path)
    seg[, CALL := derive_call(MEAN_LOG2_COPY_RATIO, gain_thr, loss_thr)]
    seg[, run_sample    := run_sample]
    seg[, patient_id    := patient_id]
    seg[, sample        := wes_id]
    seg[, wes_sample_id := wes_id]

    setcolorder(seg, c(
      "run_sample", "patient_id", "sample", "wes_sample_id",
      "CONTIG", "START", "END",
      "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO", "CALL",
      "depth", "weight", "gene"
    ))

    n_gain    <- sum(seg$CALL == "+", na.rm = TRUE)
    n_loss    <- sum(seg$CALL == "-", na.rm = TRUE)
    n_neutral <- sum(seg$CALL == "0", na.rm = TRUE)

    if (!dry_run) {
      out_parent <- dirname(out_file)
      if (!dir.exists(out_parent)) dir.create(out_parent, recursive = TRUE, showWarnings = FALSE)
      fwrite(seg, out_file)
    }
    log_msg("[Save] %s | segs=%d (gain=%d loss=%d neut=%d) | %s",
            run_sample, nrow(seg), n_gain, n_loss, n_neutral, out_file)

    summary <- rbind(summary, data.table(
      run_sample = run_sample, patient_id = patient_id, wes_sample_id = wes_id,
      cns_path = cns_path, n_segments = as.integer(nrow(seg)),
      n_gain = as.integer(n_gain), n_loss = as.integer(n_loss),
      n_neutral = as.integer(n_neutral),
      out_file = out_file, status = ifelse(dry_run, "dry_run", "saved")
    ), fill = TRUE)
  }

  if (!is.null(summary_csv) && nzchar(summary_csv)) {
    fwrite(summary, summary_csv)
    log_msg("[Summary] %s", summary_csv)
  }

  log_msg("[Done] processed=%d | saved=%d | exists=%d | missing_cns=%d",
          nrow(summary),
          sum(summary$status %in% c("saved", "dry_run")),
          sum(summary$status == "exists_skip"),
          sum(summary$status == "missing_cns"))
}

main()
