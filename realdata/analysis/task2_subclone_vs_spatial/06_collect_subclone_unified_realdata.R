#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

DEFAULT_DATASETS_ROOT <- "<DATA_ROOT>/realdata/datasets"
DEFAULT_OUT_DIR <- "<DATA_ROOT>/realdata/results/task2_subclone_vs_spatial"
DEFAULT_DATASETS <- c(
  "Ding_nature_wes",
  "metastatic_prostate_ffpe",
  "nature_prostate_ff",
  "prostate_eva_ffpe",
  "skin_tumor_ff"
)
# Per-dataset subfolder under datasets/<id>/ that contains <sample>/<run>/<tool>.
DATASET_RUN_SUBDIR <- list(
  Ding_nature_wes          = "run_links",
  metastatic_prostate_ffpe = "run_links",
  nature_prostate_ff       = "run_links",
  prostate_eva_ffpe        = "run_CNV",
  skin_tumor_ff            = "run_CNV"
)
DEFAULT_OUT_LONG <- "subclone_predictions_unified.csv"
DEFAULT_OUT_SUMMARY <- "subclone_predictions_summary.csv"
DEFAULT_PER_TOOL_OUT <- "subclone_spots_x_subclones.csv"
DEFAULT_TOOLS <- c("STARCH", "Numbat", "InferCNV", "SCEVAN", "TUSCAN", "CalicoST")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0L) stop("Cannot resolve script path.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "06_subclone_helpers_realdata.R"))

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_tool_name <- function(x) {
  k <- normalize_key(x)
  if (k == "starch") return("STARCH")
  if (k == "numbat") return("Numbat")
  if (k == "infercnv") return("InferCNV")
  if (k == "scevan") return("SCEVAN")
  if (k == "tuscan") return("TUSCAN")
  if (k == "calicost") return("CalicoST")
  NA_character_
}

log_msg <- function(...) cat(sprintf(...), "\n")

write_tool_subclone_table <- function(pred_df, tool_dir, out_name = DEFAULT_PER_TOOL_OUT) {
  if (is.null(pred_df) || nrow(pred_df) == 0L) {
    return(list(out_file = NA_character_, write_status = "no_rows"))
  }
  if (is.na(tool_dir) || !dir.exists(tool_dir)) {
    return(list(out_file = NA_character_, write_status = "tool_dir_missing"))
  }

  out_file <- file.path(tool_dir, out_name)
  dt <- as.data.table(pred_df)
  dt <- dt[, .(
    spot_id = as.character(spot_id),
    subclone = as.character(subclone)
  )]

  ok <- tryCatch(
    {
      fwrite(dt, out_file)
      TRUE
    },
    error = function(e) FALSE
  )

  list(
    out_file = out_file,
    write_status = ifelse(ok, "written", "write_failed")
  )
}

extract_for_tool <- function(tool, run_dir) {
  if (tool == "InferCNV") {
    return(extract_infercnv_with_regular_fallback(run_dir))
  }
  if (tool == "SCEVAN") {
    tool_dir <- resolve_tool_dir_by_aliases(run_dir, c("SCEVAN"))
    out <- extract_scevan_table(tool_dir)
    out$tool_dir <- tool_dir
    return(out)
  }
  if (tool == "STARCH") {
    tool_dir <- resolve_tool_dir_by_aliases(run_dir, c("STARCH"))
    out <- extract_starch_table(tool_dir)
    out$tool_dir <- tool_dir
    return(out)
  }
  if (tool == "Numbat") {
    tool_dir <- resolve_tool_dir_by_aliases(run_dir, c("Numbat"))
    out <- extract_numbat_table(tool_dir)
    out$tool_dir <- tool_dir
    return(out)
  }
  if (tool == "TUSCAN") {
    tool_dir <- resolve_tool_dir_by_aliases(run_dir, c("TUSCAN"))
    out <- extract_tuscan_table(tool_dir)
    out$tool_dir <- tool_dir
    return(out)
  }
  if (tool == "CalicoST") {
    tool_dir <- resolve_tool_dir_by_aliases(run_dir, c("CalicoST", "calicost"))
    out <- extract_calicost_table(tool_dir)
    out$tool_dir <- tool_dir
    return(out)
  }

  list(
    df = NULL,
    input_file = NA_character_,
    detail = "unsupported_tool",
    tool_dir = NA_character_
  )
}

resolve_dataset_root <- function(datasets_root, dataset) {
  subdir <- DATASET_RUN_SUBDIR[[dataset]]
  if (is.null(subdir)) {
    # Unknown dataset: try run_links first, then run_CNV.
    cand <- c(file.path(datasets_root, dataset, "run_links"),
              file.path(datasets_root, dataset, "run_CNV"))
    hit <- cand[dir.exists(cand)]
    if (length(hit) > 0L) return(hit[1])
    return(NA_character_)
  }
  path <- file.path(datasets_root, dataset, subdir)
  if (!dir.exists(path)) return(NA_character_)
  path
}

discover_runs_dataset <- function(dataset_root, samples = NULL, runs = NULL) {
  sample_dirs <- list.dirs(dataset_root, recursive = FALSE, full.names = TRUE)
  # Drop discard_* and similar housekeeping folders.
  sample_dirs <- sample_dirs[!grepl("^discard", basename(sample_dirs), ignore.case = TRUE)]
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[basename(sample_dirs) %in% samples]
  }
  if (length(sample_dirs) == 0L) {
    return(data.table(sample = character(), run = character(), run_dir = character()))
  }

  out <- vector("list", length(sample_dirs))
  oi <- 0L
  for (sdir in sample_dirs) {
    run_dirs <- list.dirs(sdir, recursive = FALSE, full.names = TRUE)
    # Only keep folders that look like run combos (i.e. contain at least one tool dir).
    if (!is.null(runs)) {
      run_dirs <- run_dirs[basename(run_dirs) %in% runs]
    }
    if (length(run_dirs) == 0L) next

    oi <- oi + 1L
    out[[oi]] <- data.table(
      sample = basename(sdir),
      run = basename(run_dirs),
      run_dir = normalizePath(run_dirs, mustWork = TRUE)
    )
  }

  if (oi == 0L) {
    return(data.table(sample = character(), run = character(), run_dir = character()))
  }
  rbindlist(out[seq_len(oi)], use.names = TRUE)
}

main <- function() {
  datasets_root <- DEFAULT_DATASETS_ROOT
  datasets <- DEFAULT_DATASETS
  out_dir <- DEFAULT_OUT_DIR
  out_long <- DEFAULT_OUT_LONG
  out_summary <- DEFAULT_OUT_SUMMARY
  samples <- NULL
  runs <- NULL
  tools <- DEFAULT_TOOLS

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L) {
    i <- 1L
    while (i <= length(args)) {
      key <- args[i]
      if (key == "--datasets_root" && i + 1L <= length(args)) {
        datasets_root <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--datasets" && i + 1L <= length(args)) {
        ds <- parse_csv_arg(args[i + 1L])
        if (length(ds) > 0L) datasets <- ds
        i <- i + 2L; next
      }
      if (key == "--out_dir" && i + 1L <= length(args)) {
        out_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_long" && i + 1L <= length(args)) {
        out_long <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_summary" && i + 1L <= length(args)) {
        out_summary <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--samples" && i + 1L <= length(args)) {
        samples <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--runs" && i + 1L <= length(args)) {
        runs <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--tools" && i + 1L <= length(args)) {
        raw_tools <- parse_csv_arg(args[i + 1L])
        mapped <- unique(vapply(raw_tools, normalize_tool_name, character(1)))
        mapped <- mapped[!is.na(mapped)]
        if (length(mapped) > 0L) tools <- mapped
        i <- i + 2L
        next
      }
      if (key %in% c("--help", "-h")) {
        cat(
          "Usage: Rscript 06_collect_subclone_unified_realdata.R [options]\n",
          "Iterates all datasets under --datasets_root and writes one merged\n",
          "subclone_predictions_unified.csv + subclone_predictions_summary.csv to --out_dir.\n",
          "Per-tool subclone_spots_x_subclones.csv files are still written next to each tool dir.\n",
          "Options:\n",
          "  --datasets_root <path>  Root containing dataset folders [default unified datasets/]\n",
          "  --datasets <csv>        Dataset whitelist (default: all 5 real datasets)\n",
          "  --out_dir <path>        Output directory for merged tables\n",
          "  --out_long <name>       Unified long CSV filename\n",
          "  --out_summary <name>    Summary CSV filename\n",
          "  --samples <csv>         Sample whitelist (applied within every dataset)\n",
          "  --runs <csv>            Run whitelist (applied within every dataset)\n",
          "  --tools <csv>           Tool list from STARCH,Numbat,InferCNV,SCEVAN,TUSCAN,CalicoST\n",
          sep = ""
        )
        quit(save = "no", status = 0L)
      }
      stop("Unknown or incomplete argument: ", key)
    }
  }

  if (!dir.exists(datasets_root)) {
    stop("datasets_root does not exist: ", datasets_root)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_long_path <- file.path(out_dir, out_long)
  out_summary_path <- file.path(out_dir, out_summary)

  log_msg("[Init] Datasets root: %s", datasets_root)
  log_msg("[Init] Datasets: %s", paste(datasets, collapse = ", "))
  log_msg("[Init] Tools: %s", paste(tools, collapse = ", "))
  log_msg("[Init] Out dir: %s", out_dir)

  long_rows <- list()
  summary_rows <- list()
  ri <- 0L

  for (ds in datasets) {
    ds_root <- resolve_dataset_root(datasets_root, ds)
    if (is.na(ds_root)) {
      log_msg("[Skip] dataset run-root not found for %s under %s", ds, datasets_root)
      next
    }

    run_tbl <- discover_runs_dataset(ds_root, samples = samples, runs = runs)
    if (nrow(run_tbl) == 0L) {
      log_msg("[Skip] no runs discovered for dataset %s (root=%s)", ds, ds_root)
      next
    }
    setorder(run_tbl, sample, run)
    log_msg("\n### Dataset: %s | runs: %d ###", ds, nrow(run_tbl))

    for (i in seq_len(nrow(run_tbl))) {
      sample_name <- run_tbl$sample[i]
      run_name <- run_tbl$run[i]
      run_dir_i <- run_tbl$run_dir[i]
      log_msg("\n=== %s / %s / %s ===", ds, sample_name, run_name)

      for (tool in tools) {
        ext <- extract_for_tool(tool, run_dir_i)
        pred <- ext$df
        tool_dir <- ext$tool_dir
        status <- if (!is.null(pred) && nrow(pred) > 0L) "ok" else ext$detail

        n_spots <- if (is.null(pred)) 0L else as.integer(nrow(pred))
        n_subclones <- if (is.null(pred)) 0L else as.integer(count_subclones(pred))
        write_out <- write_tool_subclone_table(pred, tool_dir, out_name = DEFAULT_PER_TOOL_OUT)

        if (!is.null(pred) && nrow(pred) > 0L) {
          dt <- as.data.table(pred)
          dt[, `:=`(
            dataset = ds,
            sample = sample_name,
            run = run_name,
            tool = tool,
            run_dir = run_dir_i,
            tool_dir = ifelse(is.na(tool_dir), NA_character_, tool_dir),
            input_file = ifelse(is.null(ext$input_file), NA_character_, ext$input_file)
          )]
          setcolorder(dt, c("dataset", "sample", "run", "tool", "spot_id", "subclone",
                            "run_dir", "tool_dir", "input_file"))
          long_rows[[length(long_rows) + 1L]] <- dt
        }

        ri <- ri + 1L
        summary_rows[[ri]] <- data.table(
          dataset = ds,
          sample = sample_name,
          run = run_name,
          tool = tool,
          run_dir = run_dir_i,
          tool_dir = ifelse(is.na(tool_dir), NA_character_, tool_dir),
          input_file = ifelse(is.null(ext$input_file), NA_character_, ext$input_file),
          pred_file = ifelse(is.null(write_out$out_file), NA_character_, write_out$out_file),
          n_spots = n_spots,
          n_subclones = n_subclones,
          status = status,
          detail = ext$detail,
          write_status = write_out$write_status
        )

        log_msg("  [%s] %s (spots=%d, subclones=%d, file=%s)",
                tool, status, n_spots, n_subclones, write_out$write_status)
      }
    }
  }

  if (length(summary_rows) == 0L) {
    stop("No runs processed. Check --datasets_root, --datasets and folder layout.")
  }

  summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  unified_dt <- if (length(long_rows) > 0L) {
    rbindlist(long_rows, use.names = TRUE, fill = TRUE)
  } else {
    data.table(
      dataset = character(), sample = character(), run = character(),
      tool = character(), spot_id = character(), subclone = character(),
      run_dir = character(), tool_dir = character(), input_file = character()
    )
  }

  setorder(summary_dt, dataset, sample, run, tool)
  if (nrow(unified_dt) > 0L) {
    setorder(unified_dt, dataset, sample, run, tool, spot_id)
  }

  fwrite(unified_dt, out_long_path)
  fwrite(summary_dt, out_summary_path)

  log_msg("\n[Done] Unified subclone table: %s", out_long_path)
  log_msg("[Done] Summary table: %s", out_summary_path)
  log_msg("[Done] Unified rows: %d | Summary rows: %d", nrow(unified_dt), nrow(summary_dt))
}

main()
