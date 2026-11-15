#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_spot_id <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\.", "-", x)
}

choose_col <- function(nms, candidates, required = TRUE, label = "") {
  nms_l <- tolower(nms)
  cand_l <- tolower(candidates)
  idx <- match(cand_l, nms_l)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    if (required) stop("Cannot find ", label, " column. Tried: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  nms[idx[1]]
}

map_truth <- function(label, unknown_as_na = TRUE) {
  s <- trimws(as.character(label))
  s_l <- tolower(s)

  # "possible tumor ..." should be treated as Tumor.
  is_tumor <- grepl("cancer|tumou?r|idc|malignan|possible\\s+tumou?r", s_l) | grepl("^gleason", s_l)
  is_normal <- grepl("benign|normal", s_l)
  truth_label <- ifelse(
    is_tumor, "Tumor",
    ifelse(is_normal, "Normal", ifelse(unknown_as_na, "Others", "Normal"))
  )
  truth_bin <- ifelse(truth_label == "Tumor", 1L,
                      ifelse(truth_label == "Normal", 0L, NA_integer_))
  data.table(truth_label = truth_label, truth_bin = truth_bin)
}

empty_truth <- function() {
  data.table(
    dataset_id = character(),
    sample = character(),
    spot_id = character(),
    histology_label = character(),
    truth_label = character(),
    truth_bin = integer(),
    source_file = character()
  )
}

empty_summary <- function() {
  data.table(
    dataset_id = character(),
    sample = character(),
    n_spots = integer(),
    n_tumor = integer(),
    n_normal = integer(),
    n_unknown = integer()
  )
}

main <- function() {
  dataset_root <- "<DATA_ROOT>/realdata/datasets/metastatic_prostate_ffpe"
  dataset_id <- "metastatic_prostate_ffpe"
  out_dir <- file.path(dataset_root, "evaluation", "task3_tumor_vs_histology", "truth_standardized")
  out_file <- "tumor_normal_truth.csv"
  summary_file <- "tumor_normal_truth_summary.csv"
  samples <- NULL
  overwrite <- FALSE
  unknown_as_na <- TRUE

  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--dataset_root" && i + 1L <= length(args)) {
      dataset_root <- args[i + 1L]
      i <- i + 2L; next
    }
    if (key == "--dataset_id" && i + 1L <= length(args)) {
      dataset_id <- args[i + 1L]
      i <- i + 2L; next
    }
    if (key == "--out_dir" && i + 1L <= length(args)) {
      out_dir <- args[i + 1L]
      i <- i + 2L; next
    }
    if (key == "--out_file" && i + 1L <= length(args)) {
      out_file <- args[i + 1L]
      i <- i + 2L; next
    }
    if (key == "--summary_file" && i + 1L <= length(args)) {
      summary_file <- args[i + 1L]
      i <- i + 2L; next
    }
    if (key == "--samples" && i + 1L <= length(args)) {
      samples <- parse_csv_arg(args[i + 1L])
      i <- i + 2L; next
    }
    if (key == "--overwrite") {
      overwrite <- TRUE
      i <- i + 1L; next
    }
    if (key == "--unknown_as_normal") {
      unknown_as_na <- FALSE
      i <- i + 1L; next
    }
    i <- i + 1L
  }

  out_path <- file.path(out_dir, out_file)
  sum_path <- file.path(out_dir, summary_file)
  if (file.exists(out_path) && !overwrite) {
    message("[Skip] truth file exists: ", out_path)
    return(invisible(NULL))
  }

  per_sample_dir <- file.path(dataset_root, "truth_links", "per_sample")
  sample_dirs <- if (dir.exists(per_sample_dir)) {
    d <- list.dirs(per_sample_dir, recursive = FALSE, full.names = TRUE)
    d[dir.exists(d)]
  } else {
    character()
  }
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[basename(sample_dirs) %in% samples]
  }

  out_list <- list()
  for (sd in sample_dirs) {
    sample_id <- basename(sd)
    histo_file <- file.path(sd, "histology.csv")
    if (!file.exists(histo_file)) next

    dt <- tryCatch(fread(histo_file), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0L) next

    c_spot <- choose_col(
      names(dt),
      c("Barcode", "barcode", "spot", "spot_id", "cell", "cell_id"),
      required = TRUE,
      label = "spot"
    )
    c_hist <- setdiff(names(dt), c_spot)
    if (length(c_hist) == 0L) next
    c_hist <- c_hist[1]

    x <- data.table(
      dataset_id = dataset_id,
      sample = sample_id,
      spot_id = normalize_spot_id(dt[[c_spot]]),
      histology_label = as.character(dt[[c_hist]]),
      source_file = normalizePath(histo_file, winslash = "/", mustWork = FALSE)
    )
    y <- map_truth(x$histology_label, unknown_as_na = unknown_as_na)
    x <- cbind(x, y)
    x <- x[!is.na(spot_id) & nzchar(spot_id)]

    out_list[[length(out_list) + 1L]] <- x
  }

  truth <- if (length(out_list)) {
    rbindlist(out_list, use.names = TRUE, fill = TRUE)
  } else {
    empty_truth()
  }
  setcolorder(truth, c("dataset_id", "sample", "spot_id", "histology_label", "truth_label", "truth_bin", "source_file"))

  summary_dt <- if (nrow(truth)) {
    truth[, .(
      n_spots = .N,
      n_tumor = sum(truth_bin == 1L, na.rm = TRUE),
      n_normal = sum(truth_bin == 0L, na.rm = TRUE),
      n_unknown = sum(is.na(truth_bin))
    ), by = .(dataset_id, sample)]
  } else {
    empty_summary()
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(truth, out_path)
  fwrite(summary_dt, sum_path)

  message("[Saved] ", out_path)
  message("[Saved] ", sum_path)
}

main()
