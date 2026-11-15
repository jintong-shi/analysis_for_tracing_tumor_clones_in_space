#!/usr/bin/env Rscript

# Standardize tumor/normal predictions from CNV tool outputs into a unified schema.
#
# Per tool, this script reads the canonical label/assignment file:
#   copykat : <tool_dir>/copykat_label.csv    (cols: barcode, label, is_ref)
#   scevan  : <tool_dir>/SCEVAN_label.csv     (cols: barcode, label, is_ref)
#   copyvae : <tool_dir>/clone_assignments.csv (cols: cell, clone_label, is_normal)
#   tuscan  : <tool_dir>/tumor_obj.rds.gz     (slot @tumor_cluster, named by barcode;
#                                              fallback to <run_dir>/keep_barcodes.txt
#                                              if names are missing)
#   numbat  : Numbat$new(<tool_dir>) -> nb$clone_post[, c('cell','compartment_opt')]
#             (fallback: read clone_post_2.tsv directly)
#
# Outputs (under <realdata_root>/results/task3_tumor_vs_histology/predictions_standardized/):
#   <dataset>/<sample>/<run>/<tool>/tumor_normal_standardized.csv
#   predictions_tumor_normal_long.csv
#   predictions_processing_summary.csv

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

normalize_spot_id <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\.", "-", x)
}

norm_tool <- function(x) {
  y <- tolower(trimws(as.character(x)))
  if (y %in% c("copykat")) return("copykat")
  if (y %in% c("scevan"))  return("scevan")
  if (y %in% c("copyvae")) return("copyvae")
  if (y %in% c("tuscan"))  return("tuscan")
  if (y %in% c("numbat"))  return("numbat")
  NA_character_
}

label_to_bin <- function(label) {
  s <- tolower(trimws(as.character(label)))
  bin <- rep(NA_integer_, length(s))
  # Substring match — copykat sometimes emits compound labels like
  # "c1:diploid:low.conf" / "c2:aneuploid:low.conf".
  is_tumor  <- grepl("aneuploid|tumou?r|malignan|cancer", s)
  is_normal <- grepl("diploid|^normal$|benign", s)
  is_unknown <- grepl("not[._ ]?defined|unknown|unassigned|na", s) & !is_tumor & !is_normal
  bin[is_tumor & !is_normal] <- 1L
  bin[is_normal & !is_tumor] <- 0L
  out_label <- ifelse(bin == 1L, "Tumor",
                      ifelse(bin == 0L, "Normal", NA_character_))
  data.table(pred_label = out_label, pred_bin = bin)
}

# ---- per-tool extractors ----------------------------------------------------

extract_copykat <- function(tool_dir) {
  f <- file.path(tool_dir, "copykat_label.csv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  if (!all(c("barcode", "label") %in% names(dt))) return(NULL)
  if ("is_ref" %in% names(dt)) dt <- dt[!isTRUE(is_ref) & is_ref != TRUE]
  if (nrow(dt) == 0L) return(NULL)
  z <- label_to_bin(dt$label)
  data.table(
    spot_id = normalize_spot_id(dt$barcode),
    z,
    source_type = "copykat_label_csv"
  )
}

extract_scevan <- function(tool_dir) {
  f <- file.path(tool_dir, "SCEVAN_label.csv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  if (!all(c("barcode", "label") %in% names(dt))) return(NULL)
  if ("is_ref" %in% names(dt)) dt <- dt[!isTRUE(is_ref) & is_ref != TRUE]
  if (nrow(dt) == 0L) return(NULL)
  z <- label_to_bin(dt$label)
  data.table(
    spot_id = normalize_spot_id(dt$barcode),
    z,
    source_type = "scevan_label_csv"
  )
}

extract_copyvae <- function(tool_dir) {
  f <- file.path(tool_dir, "clone_assignments.csv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  if (!all(c("cell", "is_normal") %in% names(dt))) return(NULL)

  is_norm_chr <- tolower(trimws(as.character(dt$is_normal)))
  bin <- ifelse(is_norm_chr %in% c("true", "1", "t"), 0L,
                ifelse(is_norm_chr %in% c("false", "0", "f"), 1L, NA_integer_))
  lab <- ifelse(bin == 1L, "Tumor", ifelse(bin == 0L, "Normal", NA_character_))
  data.table(
    spot_id = normalize_spot_id(dt$cell),
    pred_label = lab,
    pred_bin = bin,
    source_type = "copyvae_clone_assignments"
  )
}

extract_tuscan <- function(tool_dir, run_dir) {
  f <- file.path(tool_dir, "tumor_obj.rds.gz")
  if (!file.exists(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  labs <- tryCatch(obj@tumor_cluster, error = function(e) NULL)
  if (is.null(labs) || length(labs) == 0L) return(NULL)

  cells <- names(labs)
  if (is.null(cells) || any(!nzchar(cells))) {
    img <- tryCatch(obj@image, error = function(e) NULL)
    if (!is.null(img) && is.data.frame(img) && "barcode" %in% names(img)) {
      cand <- as.character(img$barcode)
      if (length(cand) == length(labs)) cells <- cand
    }
    if (is.null(cells) || length(cells) != length(labs)) {
      cnv_cols <- tryCatch(colnames(obj@cnv.data), error = function(e) NULL)
      if (!is.null(cnv_cols) && length(cnv_cols) == length(labs)) cells <- as.character(cnv_cols)
    }
    # Fallback: order matches obj@annotation, otherwise keep_barcodes.txt.
    anno <- tryCatch(obj@annotation, error = function(e) NULL)
    if ((is.null(cells) || length(cells) != length(labs)) && !is.null(anno)) {
      cand <- if (!is.null(rownames(anno))) rownames(anno) else
              if (!is.null(names(anno))) names(anno) else
              if (is.data.frame(anno) && "barcode" %in% names(anno)) anno$barcode else NULL
      if (!is.null(cand) && length(cand) == length(labs)) cells <- as.character(cand)
    }
    if (is.null(cells) || length(cells) != length(labs)) {
      kb <- file.path(run_dir, "keep_barcodes.txt")
      if (file.exists(kb)) {
        bc <- tryCatch(readLines(kb, warn = FALSE), error = function(e) character())
        bc <- bc[nzchar(bc)]
        if (length(bc) == length(labs)) cells <- bc
      }
    }
  }
  if (is.null(cells) || length(cells) != length(labs)) return(NULL)

  z <- label_to_bin(as.character(labs))
  data.table(
    spot_id = normalize_spot_id(cells),
    z,
    source_type = "tuscan_tumor_obj_rds"
  )
}

extract_numbat <- function(tool_dir) {
  # numbat sometimes produces no CNV-positive iterations and writes a sentinel.
  # In that case clone_post / joint_post / etc. are absent and there is nothing
  # to standardize.
  if (file.exists(file.path(tool_dir, "NUMBAT.NO_CNV"))) return(NULL)

  cp <- NULL
  src <- NA_character_

  # Preferred path: use the numbat package.
  # Numbat$new() greedily fread()s a long list of optional files and prints
  # benign warnings for any that are missing — silence them.
  if (requireNamespace("numbat", quietly = TRUE)) {
    nb <- tryCatch(
      suppressWarnings(suppressMessages(numbat::Numbat$new(tool_dir))),
      error = function(e) NULL
    )
    if (!is.null(nb)) {
      cp <- tryCatch(as.data.frame(nb$clone_post), error = function(e) NULL)
      if (!is.null(cp)) src <- "numbat_clone_post"
    }
  }

  # Fallback: read the latest clone_post_<i>.tsv file.
  if (is.null(cp) || !all(c("cell", "compartment_opt") %in% names(cp))) {
    cps <- list.files(tool_dir, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
    if (length(cps) == 0L) return(NULL)
    idx <- as.integer(sub("^clone_post_(\\d+)\\.tsv$", "\\1", basename(cps)))
    cp <- tryCatch(fread(cps[which.max(idx)]), error = function(e) NULL)
    if (is.null(cp)) return(NULL)
    src <- "numbat_clone_post_tsv"
  }
  if (!all(c("cell", "compartment_opt") %in% names(cp))) return(NULL)

  z <- label_to_bin(cp$compartment_opt)
  data.table(
    spot_id = normalize_spot_id(cp$cell),
    z,
    source_type = src
  )
}

extract_tool <- function(tool, tool_dir, run_dir) {
  switch(tool,
         copykat = extract_copykat(tool_dir),
         scevan  = extract_scevan(tool_dir),
         copyvae = extract_copyvae(tool_dir),
         tuscan  = extract_tuscan(tool_dir, run_dir),
         numbat  = extract_numbat(tool_dir),
         NULL)
}

# ---- main -------------------------------------------------------------------

dataset_id_from_root <- function(path) {
  basename(normalizePath(path, winslash = "/", mustWork = FALSE))
}

main <- function() {
  dataset_root <- NULL
  realdata_root <- NULL
  dataset_ids <- NULL
  truth_csv <- NULL
  out_dir <- NULL
  tools <- c("copykat", "scevan", "copyvae", "tuscan", "numbat")
  samples <- NULL
  runs <- NULL
  sample_overview <- NULL

  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--dataset_root" && i + 1L <= length(args)) {
      dataset_root <- args[i + 1L]
      dataset_ids <- unique(c(dataset_ids, dataset_id_from_root(dataset_root)))
      i <- i + 2L; next
    }
    if (key == "--realdata_root" && i + 1L <= length(args)) {
      realdata_root <- args[i + 1L]; i <- i + 2L; next
    }
    if (key %in% c("--dataset_id", "--dataset_ids", "--datasets") && i + 1L <= length(args)) {
      dataset_ids <- unique(c(dataset_ids, parse_csv_arg(args[i + 1L]))); i <- i + 2L; next
    }
    if (key == "--truth_csv" && i + 1L <= length(args)) {
      truth_csv <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--out_dir" && i + 1L <= length(args)) {
      out_dir <- args[i + 1L]; i <- i + 2L; next
    }
    if (key == "--tools" && i + 1L <= length(args)) {
      tools <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
    }
    if (key == "--samples" && i + 1L <= length(args)) {
      samples <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
    }
    if (key == "--runs" && i + 1L <= length(args)) {
      runs <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
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
  if (is.null(out_dir)) out_dir <- task3_results_dir(realdata_root)

  task_samples <- read_task3_sample_table(
    sample_overview = sample_overview,
    realdata_root = realdata_root,
    dataset_ids = dataset_ids,
    samples = samples
  )
  if (!nrow(task_samples)) stop("No Task3 histology-available samples selected.")

  tool_set <- unique(na.omit(vapply(tools, norm_tool, character(1))))
  if (length(tool_set) == 0L) stop("No valid tools requested.")

  pred_dir <- file.path(out_dir, "predictions_standardized")
  dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  sm <- list()

  for (ds_id in unique(task_samples$dataset_id)) {
    dataset_id <- ds_id
    dataset_id_value <- ds_id
    ds_samples <- task_samples[dataset_id == ds_id, sample]
    ds_root <- if (!is.null(dataset_root) && dataset_id_from_root(dataset_root) == dataset_id) {
      dataset_root
    } else {
      task3_dataset_root(realdata_root, dataset_id)
    }
    if (!dir.exists(ds_root)) stop("Missing dataset root for ", dataset_id, ": ", ds_root)

    truth_csv_ds <- if (!is.null(truth_csv) && length(unique(task_samples$dataset_id)) == 1L) {
      truth_csv
    } else {
      task3_truth_csv(ds_root)
    }
    run_root <- task3_run_root(ds_root)

    truth <- read_task3_truth_for_dataset(
      truth_csv = truth_csv_ds,
      task_samples = task_samples,
      dataset_id = dataset_id
    )
    truth[, dataset_id := dataset_id_value]
    truth[, sample := normalize_sample_id(sample)]
    truth[, spot_id := normalize_spot_id(spot_id)]
    truth <- truth[sample %in% ds_samples]

    message("[Sample overview] ", dataset_id, ": selected ", length(ds_samples),
            " histology-available samples")

    for (sample_id in ds_samples) {
      sample_dir <- file.path(run_root, sample_id)
      n_truth_sample <- truth[sample == sample_id, .N]
      if (!dir.exists(sample_dir)) {
        message("[Skip] ", dataset_id, "/", sample_id, ": no run directory")
        sm[[length(sm) + 1L]] <- data.table(
          dataset_id = dataset_id, sample = sample_id, run = NA_character_, tool = NA_character_,
          n_truth = n_truth_sample, n_pred = 0L, n_used = 0L,
          status = "sample_dir_missing", tool_dir = NA_character_
        )
        next
      }

      run_dirs <- list.files(sample_dir, full.names = TRUE)
      run_dirs <- run_dirs[dir.exists(run_dirs)]
      if (!is.null(runs)) run_dirs <- run_dirs[basename(run_dirs) %in% runs]

      if (n_truth_sample == 0L) {
        message("[Skip] ", dataset_id, "/", sample_id, ": no truth rows")
        sm[[length(sm) + 1L]] <- data.table(
          dataset_id = dataset_id, sample = sample_id, run = NA_character_, tool = NA_character_,
          n_truth = 0L, n_pred = 0L, n_used = 0L,
          status = "truth_missing", tool_dir = NA_character_
        )
        next
      }

      for (rlnk in run_dirs) {
        run_id <- basename(rlnk)
        run_dir <- normalizePath(rlnk, winslash = "/", mustWork = FALSE)
        if (!dir.exists(run_dir)) next

        tool_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
        tool_dirs <- tool_dirs[dir.exists(tool_dirs)]

        for (td in tool_dirs) {
          tool_id <- norm_tool(basename(td))
          if (is.na(tool_id) || !(tool_id %in% tool_set)) next

          pred <- tryCatch(extract_tool(tool_id, td, run_dir),
                           error = function(e) {
                             message("[Error] ", dataset_id, "/", sample_id, "/",
                                     run_id, "/", tool_id, ": ", conditionMessage(e))
                             NULL
                           })

          if (is.null(pred) || nrow(pred) == 0L) {
            sm[[length(sm) + 1L]] <- data.table(
              dataset_id = dataset_id, sample = sample_id, run = run_id, tool = tool_id,
              n_truth = n_truth_sample, n_pred = 0L, n_used = 0L,
              status = "prediction_missing", tool_dir = td
            )
            next
          }

          pred[, spot_id := normalize_spot_id(spot_id)]
          pred <- unique(pred[!is.na(spot_id) & nzchar(spot_id),
                              .(spot_id, pred_label, pred_bin, source_type)])

          tt <- truth[sample == sample_id,
                      .(dataset_id, sample, spot_id, histology_label, truth_label, truth_bin)]
          mg <- merge(tt, pred, by = "spot_id", all = FALSE)
          mg[, `:=`(
            dataset_id = dataset_id_value,
            sample = sample_id,
            run = run_id,
            tool = tool_id,
            tool_dir = td
          )]

          out_subdir <- file.path(pred_dir, dataset_id, sample_id, run_id, tool_id)
          dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
          fwrite(mg, file.path(out_subdir, "tumor_normal_standardized.csv"))
          rows[[length(rows) + 1L]] <- mg

          sm[[length(sm) + 1L]] <- data.table(
            dataset_id = dataset_id, sample = sample_id, run = run_id, tool = tool_id,
            n_truth = nrow(tt), n_pred = nrow(pred),
            n_used = sum(!is.na(mg$truth_bin) & !is.na(mg$pred_bin)),
            status = "ok", tool_dir = td
          )
        }
      }
    }
  }

  all_pred <- if (length(rows)) rbindlist(rows, use.names = TRUE, fill = TRUE) else
    data.table(spot_id = character(), dataset_id = character(), sample = character(),
               histology_label = character(), truth_label = character(), truth_bin = integer(),
               pred_label = character(), pred_bin = integer(), source_type = character(),
               run = character(), tool = character(), tool_dir = character())

  all_sm <- if (length(sm)) rbindlist(sm, use.names = TRUE, fill = TRUE) else
    data.table(dataset_id = character(), sample = character(), run = character(), tool = character(),
               n_truth = integer(), n_pred = integer(), n_used = integer(),
               status = character(), tool_dir = character())

  out_long <- file.path(pred_dir, "predictions_tumor_normal_long.csv")
  out_summary <- file.path(pred_dir, "predictions_processing_summary.csv")
  fwrite(all_pred, out_long)
  fwrite(all_sm, out_summary)

  message("[Saved] ", out_long)
  message("[Saved] ", out_summary)
}

main()
