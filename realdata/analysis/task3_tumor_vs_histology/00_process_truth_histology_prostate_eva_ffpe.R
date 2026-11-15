#!/usr/bin/env Rscript

# Step-0 truth builder for the `prostate_eva_ffpe` dataset.
#
# Histology truth lives in `truth_links/histology/<file>_three_categories.csv`.
# Several files cover one section per Visium sample; two combined-slide files
# (PCa20153_C1_20128_C1, PCa20130_C1_20272_C2) stitch two sub-sections onto a
# single Visium capture, so we keep both sub-sections as truth for the matching
# `P10_P3` / `P4_P13` runs.
#
# As a convenience this script also creates a `run_links` -> `run_CNV` symlink
# and per-sample `truth_links/per_sample/<sample>/spatial/tissue_positions.csv`
# symlinks so that scripts/01..04 (which expect the nature_prostate_ff layout)
# can run unchanged.

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

# Sample -> (histology file root, vector of sub-section ids to keep)
# Mapping derived from public_data/.../section_name_mapping.xlsx (sheet1).
# `sections = NULL` means: file has only one section, keep all rows.
sample_truth_map <- list(
  P25     = list(file = "A1_20033_three_categories.csv",                sections = NULL),
  P11_1   = list(file = "A1_20216_C1_three_categories.csv",             sections = NULL),
  P15_2   = list(file = "B1_19617_C2_three_categories.csv",             sections = NULL),
  P2_2    = list(file = "B1_20111_C2_three_categories.csv",             sections = NULL),
  P4_2    = list(file = "D1_20130_C2_three_categories.csv",             sections = NULL),
  P10_2   = list(file = "D1_20153_C2_three_categories.csv",             sections = NULL),
  P10_P3  = list(file = "PCa20153_C1_20128_C1_three_categories.csv",    sections = c("20153-1", "20128-1")),
  P4_P13  = list(file = "PCa20130_C1_20272_C2_three_categories.csv",    sections = c("20272-2", "20130-1"))
)
# The map keys use underscores; convert back to the canonical run_CNV folder
# names ("P11_1" -> "P11-1") via this lookup.
sample_folder_for <- c(
  P25 = "P25", P11_1 = "P11-1", P15_2 = "P15-2", P2_2 = "P2-2",
  P4_2 = "P4-2", P10_2 = "P10-2", P10_P3 = "P10_P3", P4_P13 = "P4_P13"
)

map_truth_from_final <- function(final) {
  s <- trimws(as.character(final))
  s_l <- tolower(s)
  truth_label <- ifelse(s_l == "tumor", "Tumor",
                 ifelse(s_l == "normal", "Normal", "Others"))
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

ensure_symlink <- function(target, link_path) {
  if (file.exists(link_path) || file.symlink(target, link_path)) {
    return(invisible(TRUE))
  }
  message("[Warn] could not create symlink: ", link_path, " -> ", target)
  invisible(FALSE)
}

setup_run_and_spatial_links <- function(dataset_root, samples_done) {
  run_links <- file.path(dataset_root, "run_links")
  run_cnv <- file.path(dataset_root, "run_CNV")
  if (!file.exists(run_links) && dir.exists(run_cnv)) {
    if (file.symlink(run_cnv, run_links)) {
      message("[Linked] ", run_links, " -> ", run_cnv)
    }
  }

  per_sample <- file.path(dataset_root, "truth_links", "per_sample")
  dir.create(per_sample, recursive = TRUE, showWarnings = FALSE)
  for (s in samples_done) {
    src <- file.path(dataset_root, "source_links", paste0("spaceranger_", s),
                     "outs", "spatial", "tissue_positions.csv")
    if (!file.exists(src)) {
      message("[Warn] spaceranger spatial missing for ", s, ": ", src)
      next
    }
    dst_dir <- file.path(per_sample, s, "spatial")
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    dst <- file.path(dst_dir, "tissue_positions.csv")
    if (!file.exists(dst)) {
      if (file.symlink(src, dst)) {
        message("[Linked] ", dst, " -> ", src)
      } else {
        message("[Warn] could not create symlink: ", dst)
      }
    }
  }
}

main <- function() {
  dataset_root <- "<DATA_ROOT>/realdata/datasets/prostate_eva_ffpe"
  dataset_id <- "prostate_eva_ffpe"
  out_dir <- file.path(dataset_root, "evaluation", "task3_tumor_vs_histology", "truth_standardized")
  out_file <- "tumor_normal_truth.csv"
  summary_file <- "tumor_normal_truth_summary.csv"
  mapping_file <- "label_mapping_summary.csv"
  samples <- NULL
  overwrite <- FALSE
  do_links <- TRUE

  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--dataset_root" && i + 1L <= length(args)) {
      dataset_root <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--dataset_id" && i + 1L <= length(args)) {
      dataset_id <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--out_dir" && i + 1L <= length(args)) {
      out_dir <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--out_file" && i + 1L <= length(args)) {
      out_file <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--summary_file" && i + 1L <= length(args)) {
      summary_file <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--mapping_file" && i + 1L <= length(args)) {
      mapping_file <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--samples" && i + 1L <= length(args)) {
      samples <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
    }
    if (key == "--overwrite") {
      overwrite <- TRUE; i <- i + 1L; next
    }
    if (key == "--no_links") {
      do_links <- FALSE; i <- i + 1L; next
    }
    i <- i + 1L
  }

  out_path <- file.path(out_dir, out_file)
  sum_path <- file.path(out_dir, summary_file)
  map_path <- file.path(out_dir, mapping_file)
  histo_dir <- file.path(dataset_root, "truth_links", "histology")

  if (file.exists(out_path) && !overwrite) {
    message("[Skip] truth file exists: ", out_path)
    if (do_links) {
      setup_run_and_spatial_links(dataset_root, unname(sample_folder_for))
    }
    return(invisible(NULL))
  }

  out_list <- list()
  for (key in names(sample_truth_map)) {
    sample_id <- sample_folder_for[[key]]
    if (!is.null(samples) && !(sample_id %in% samples)) next

    cfg <- sample_truth_map[[key]]
    histo_file <- file.path(histo_dir, cfg$file)
    if (!file.exists(histo_file)) {
      message("[Warn] missing truth file for ", sample_id, ": ", histo_file)
      next
    }

    dt <- tryCatch(fread(histo_file), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0L) next

    c_spot <- choose_col(names(dt),
                         c("Barcode", "barcode", "spot", "spot_id"),
                         required = TRUE, label = "spot")
    c_hist <- choose_col(names(dt),
                         c("Histology", "histology", "histology_label"),
                         required = TRUE, label = "histology")
    c_final <- choose_col(names(dt),
                          c("Final_Annotations", "final_annotation", "final_annotations"),
                          required = TRUE, label = "Final_Annotations")
    c_sub <- choose_col(names(dt), c("sample_id"), required = FALSE)

    if (!is.null(cfg$sections) && (is.na(c_sub) || !nzchar(c_sub))) {
      stop("Truth file for ", sample_id, " is expected to have a `sample_id` ",
           "column for sub-section filtering: ", histo_file)
    }

    keep <- rep(TRUE, nrow(dt))
    if (!is.null(cfg$sections) && !is.na(c_sub) && nzchar(c_sub)) {
      keep <- as.character(dt[[c_sub]]) %in% cfg$sections
    }

    dt2 <- dt[keep]
    if (nrow(dt2) == 0L) {
      message("[Warn] no rows after sub-section filter for ", sample_id)
      next
    }

    x <- data.table(
      dataset_id = dataset_id,
      sample = sample_id,
      spot_id = normalize_spot_id(dt2[[c_spot]]),
      histology_label = as.character(dt2[[c_hist]]),
      source_file = normalizePath(histo_file, winslash = "/", mustWork = FALSE)
    )
    y <- map_truth_from_final(dt2[[c_final]])
    x <- cbind(x, y)
    x <- x[!is.na(spot_id) & nzchar(spot_id)]

    out_list[[length(out_list) + 1L]] <- x
  }

  truth <- if (length(out_list)) {
    rbindlist(out_list, use.names = TRUE, fill = TRUE)
  } else {
    empty_truth()
  }
  setcolorder(truth, c("dataset_id", "sample", "spot_id", "histology_label",
                       "truth_label", "truth_bin", "source_file"))

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

  mapping_dt <- if (nrow(truth)) {
    truth[, .N, by = .(histology_label, truth_label)][order(truth_label, -N, histology_label)]
  } else {
    data.table(histology_label = character(), truth_label = character(), N = integer())
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(truth, out_path)
  fwrite(summary_dt, sum_path)
  fwrite(mapping_dt, map_path)

  message("[Saved] ", out_path)
  message("[Saved] ", sum_path)
  message("[Saved] ", map_path)

  if (do_links) {
    samples_done <- unique(truth$sample)
    setup_run_and_spatial_links(dataset_root, samples_done)
  }
}

main()
