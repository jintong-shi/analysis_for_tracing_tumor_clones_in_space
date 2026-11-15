#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  getwd()
}

DEFAULT_ROOT_DIR <- NULL
DEFAULT_DATASETS_DIR <- normalizePath(file.path(get_script_dir(), "..", "..", "..", "datasets"), mustWork = FALSE)
DEFAULT_DECON_RUN <- "Decon_no-ref"
DEFAULT_DECON_FALLBACK_RUN <- "Decon_StromaImm-ref"
DEFAULT_DECON_KEEP_FILE <- "keep_barcodes.txt"
DEFAULT_INPUT_BASENAMES <- c("processed_log2_spots_x_genes.csv", "processed_state_spots_x_genes.csv")

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

dataset_has_task1_truth <- function(dataset_dir) {
  truth_dirs <- file.path(dataset_dir, "truth_links", c("wes_cnv", "wes_cnv_tumor_only"))
  any(dir.exists(truth_dirs))
}

discover_dataset_roots <- function(datasets_dir, datasets = NULL) {
  if (!dir.exists(datasets_dir)) stop("datasets_dir does not exist: ", datasets_dir)

  dataset_dirs <- list.dirs(datasets_dir, recursive = FALSE, full.names = TRUE)
  dataset_dirs <- dataset_dirs[file.info(dataset_dirs)$isdir]
  if (!is.null(datasets)) {
    dataset_dirs <- dataset_dirs[normalize_key(basename(dataset_dirs)) %in% normalize_key(datasets)]
  }

  rows <- list()
  for (dataset_dir in dataset_dirs) {
    run_candidates <- c(
      file.path(dataset_dir, "run_links", "run_CNV"),
      file.path(dataset_dir, "run_links"),
      file.path(dataset_dir, "run_CNV")
    )
    run_candidates <- run_candidates[dir.exists(run_candidates)]
    if (length(run_candidates) == 0L) next
    if (is.null(datasets) && !dataset_has_task1_truth(dataset_dir)) next

    rows[[length(rows) + 1L]] <- data.table(
      dataset = basename(dataset_dir),
      root_dir = normalizePath(run_candidates[1], mustWork = TRUE)
    )
  }

  if (length(rows) == 0L) return(data.table(dataset = character(0), root_dir = character(0)))
  rbindlist(rows, fill = TRUE)
}

dataset_name_from_root <- function(root_dir) {
  root_dir <- normalizePath(root_dir, mustWork = FALSE)
  base <- basename(root_dir)
  parent <- basename(dirname(root_dir))
  if (base == "run_CNV" && parent == "run_links") return(basename(dirname(dirname(root_dir))))
  if (base %in% c("run_CNV", "run_links")) return(basename(dirname(root_dir)))
  basename(dirname(root_dir))
}

normalize_spot_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  x
}

log_msg <- function(...) cat(sprintf(...), "\n")

choose_col <- function(nms, candidates, required = FALSE, label = "") {
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

read_spot_vector <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  x <- normalize_spot_ids(x)
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

read_tumor_spots_from_histo <- function(histo_dir) {
  # Preferred explicit files for future histology integration.
  candidate_spot_files <- c(
    "tumor_spots.txt",
    "tumor_barcodes.txt",
    "tumor_spot_ids.txt",
    "tumour_spots.txt",
    "keep_barcodes_tumor.txt",
    "keep_tumor_barcodes.txt",
    "keep_barcodes.txt"
  )
  for (fn in candidate_spot_files) {
    f <- file.path(histo_dir, fn)
    if (file.exists(f)) {
      spots <- read_spot_vector(f)
      if (length(spots) > 0L) return(spots)
    }
  }

  # Flexible parser for table annotations (csv/tsv/txt).
  tbl_files <- list.files(
    histo_dir,
    pattern = "\\.(csv|tsv|txt)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(tbl_files) == 0L) return(character(0))

  for (f in tbl_files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0L) next

    nms <- names(dt)
    spot_col <- choose_col(
      nms,
      c("barcode", "barcodes", "spot", "spot_id", "spotid", "cell", "cell_id"),
      required = FALSE
    )
    if (is.na(spot_col)) next

    label_col <- choose_col(
      nms,
      c("label", "class", "annotation", "histology", "celltype", "region"),
      required = FALSE
    )
    flag_col <- choose_col(
      nms,
      c("tumor", "is_tumor", "is_tumour", "epi", "is_epi", "epithelial"),
      required = FALSE
    )

    keep <- rep(FALSE, nrow(dt))
    if (!is.na(label_col)) {
      lab <- normalize_key(dt[[label_col]])
      keep <- keep | grepl("tumou?r|malignan|epi|epithelial", lab)
    }
    if (!is.na(flag_col)) {
      fv <- tolower(trimws(as.character(dt[[flag_col]])))
      keep <- keep | fv %in% c("1", "true", "t", "yes", "y", "tumor", "tumour", "epi", "epithelial")
    }

    if (!any(keep)) next
    spots <- normalize_spot_ids(dt[[spot_col]][keep])
    spots <- spots[!is.na(spots) & nzchar(spots)]
    if (length(spots) > 0L) return(unique(spots))
  }

  character(0)
}

detect_spot_selection <- function(sample_dir, decon_run, decon_keep_file, histo_fallback_to_decon = FALSE) {
  subdirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
  subdirs <- subdirs[file.info(subdirs)$isdir]
  histo_dirs <- subdirs[grepl("^Histo", basename(subdirs), ignore.case = TRUE)]

  if (length(histo_dirs) > 0L) {
    histo_dir <- sort(histo_dirs)[1]
    spots <- read_tumor_spots_from_histo(histo_dir)
    if (length(spots) > 0L) {
      return(list(
        source = paste0("histology:", basename(histo_dir)),
        spots = spots
      ))
    }

    if (!histo_fallback_to_decon) {
      return(list(
        source = paste0("histology_found_no_tumor_spots:", basename(histo_dir)),
        spots = character(0)
      ))
    }
  }

  decon_keep <- file.path(sample_dir, decon_run, decon_keep_file)
  spots <- read_spot_vector(decon_keep)
  if (length(spots) > 0L) {
    return(list(
      source = paste0("decon_keep_barcodes:", decon_run),
      spots = spots
    ))
  }

  if (!identical(decon_run, DEFAULT_DECON_FALLBACK_RUN)) {
    fallback_keep <- file.path(sample_dir, DEFAULT_DECON_FALLBACK_RUN, decon_keep_file)
    fallback_spots <- read_spot_vector(fallback_keep)
    if (length(fallback_spots) > 0L) {
      return(list(
        source = paste0("decon_keep_barcodes:", DEFAULT_DECON_FALLBACK_RUN),
        spots = fallback_spots
      ))
    }
  }

  list(
    source = paste0("decon_keep_barcodes:", decon_run),
    spots = spots
  )
}

read_processed_matrix <- function(path) {
  dt <- tryCatch(fread(path, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L || ncol(dt) < 2L) return(NULL)

  rn <- normalize_spot_ids(dt[[1]])
  keep <- !is.na(rn) & nzchar(rn)
  if (!any(keep)) return(NULL)

  rn <- rn[keep]
  mat <- as.matrix(dt[keep, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- rn
  colnames(mat) <- colnames(dt)[-1]

  if (anyDuplicated(rownames(mat)) > 0L) {
    sums <- rowsum(mat, group = rownames(mat), reorder = FALSE, na.rm = TRUE)
    cnt <- rowsum(!is.na(mat), group = rownames(mat), reorder = FALSE)
    mat <- sums / pmax(cnt, 1)
  }

  mat
}

build_pseudobulk_gene_x_1 <- function(spots_x_genes, selected_spots) {
  selected_spots <- unique(normalize_spot_ids(selected_spots))
  selected_spots <- selected_spots[!is.na(selected_spots) & nzchar(selected_spots)]
  if (length(selected_spots) == 0L) {
    return(list(gene_x_1 = NULL, matched_spots = character(0)))
  }

  idx <- match(selected_spots, rownames(spots_x_genes), nomatch = 0L)
  idx <- unique(idx[idx > 0L])
  if (length(idx) == 0L) {
    return(list(gene_x_1 = NULL, matched_spots = character(0)))
  }

  sub <- spots_x_genes[idx, , drop = FALSE]
  n_non_na <- colSums(!is.na(sub))
  sums <- colSums(sub, na.rm = TRUE)
  vec <- sums / pmax(n_non_na, 1L)
  vec[n_non_na == 0L] <- NA_real_

  out <- matrix(vec, ncol = 1L, dimnames = list(colnames(spots_x_genes), "pseudobulk"))
  list(gene_x_1 = out, matched_spots = rownames(sub))
}

make_output_basename <- function(in_basename) {
  out <- sub("^processed_", "pseudobulk_tumor_epi_", in_basename)
  out <- sub("_spots_x_genes", "_gene_x_1", out)
  out
}

# Map a tool folder to its evaluated methods. infercnv and numbat each expose two
# methods (expr + cnv) from a single folder; every other tool is a single method
# keyed on the folder name, reading the default input basenames.
tool_method_specs <- function(tool_name, input_basenames) {
  key <- normalize_key(tool_name)
  if (key == "infercnv") {
    return(list(
      list(tool = "infercnv_expr", in_basename = "processed_infercnv_expr_spots_x_genes.csv"),
      list(tool = "infercnv_cnv",  in_basename = "processed_infercnv_cnv_spots_x_genes.csv")
    ))
  }
  if (key == "numbat") {
    return(list(
      list(tool = "numbat_expr", in_basename = "processed_numbat_expr_spots_x_genes.csv"),
      list(tool = "numbat_cnv",  in_basename = "processed_numbat_cnv_spots_x_genes.csv")
    ))
  }
  lapply(input_basenames, function(b) list(tool = tool_name, in_basename = b))
}

process_one_tool <- function(sample_name, run_name, method_specs, tool_dir, selected_spots, spot_source,
                             overwrite = FALSE, dry_run = FALSE) {
  rows <- list()

  for (spec in method_specs) {
    in_basename <- spec$in_basename
    tool_name <- spec$tool
    in_file <- file.path(tool_dir, in_basename)
    out_basename <- make_output_basename(in_basename)
    out_file <- file.path(tool_dir, out_basename)

    row <- data.table(
      sample = sample_name,
      run = run_name,
      tool = tool_name,
      input_file = in_file,
      output_file = out_file,
      spot_source = spot_source,
      n_selected_spots = as.integer(length(selected_spots)),
      n_spots_in_matrix = NA_integer_,
      n_matched_spots = NA_integer_,
      n_genes = NA_integer_,
      status = "init"
    )

    if (!file.exists(in_file)) {
      row$status <- "input_missing"
      rows[[length(rows) + 1L]] <- row
      next
    }

    if (file.exists(out_file) && !overwrite) {
      row$status <- "exists_skip"
      rows[[length(rows) + 1L]] <- row
      next
    }

    mat <- read_processed_matrix(in_file)
    if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      row$status <- "matrix_parse_failed"
      rows[[length(rows) + 1L]] <- row
      next
    }
    row$n_spots_in_matrix <- as.integer(nrow(mat))
    row$n_genes <- as.integer(ncol(mat))

    pb <- build_pseudobulk_gene_x_1(mat, selected_spots)
    row$n_matched_spots <- as.integer(length(pb$matched_spots))
    if (is.null(pb$gene_x_1)) {
      row$status <- "no_matched_spots"
      rows[[length(rows) + 1L]] <- row
      next
    }

    if (!dry_run) {
      write.csv(pb$gene_x_1, out_file, quote = FALSE)
    }
    row$status <- if (dry_run) "dry_run" else "saved"
    rows[[length(rows) + 1L]] <- row
  }

  rbindlist(rows, fill = TRUE)
}

process_root <- function(dataset_name, root_dir, samples, strategies, tools, decon_run, decon_keep_file,
                         input_basenames, overwrite, dry_run, histo_fallback_to_decon) {
  if (!dir.exists(root_dir)) stop("root_dir does not exist: ", root_dir)
  if (length(input_basenames) == 0L) stop("No input basenames provided.")

  sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[file.info(sample_dirs)$isdir]
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[normalize_key(basename(sample_dirs)) %in% normalize_key(samples)]
  }

  log_msg("[Init] dataset = %s | root_dir = %s", dataset_name, root_dir)
  log_msg("[Init] sample_count = %d", length(sample_dirs))
  log_msg("[Init] decon_run = %s | decon_keep_file = %s", decon_run, decon_keep_file)
  log_msg("[Init] input_basenames = %s", paste(input_basenames, collapse = ", "))
  log_msg("[Init] overwrite = %s | dry_run = %s", as.character(overwrite), as.character(dry_run))
  log_msg("[Init] histo_fallback_to_decon = %s", as.character(histo_fallback_to_decon))

  summary_rows <- list()
  row_i <- 1L

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    sel <- detect_spot_selection(
      sample_dir = sample_dir,
      decon_run = decon_run,
      decon_keep_file = decon_keep_file,
      histo_fallback_to_decon = histo_fallback_to_decon
    )
    selected_spots <- sel$spots
    spot_source <- sel$source

    log_msg("[Sample] %s | spot_source=%s | n_selected=%d", sample_name, spot_source, length(selected_spots))

    run_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[file.info(run_dirs)$isdir]
    if (!is.null(strategies)) {
      run_dirs <- run_dirs[normalize_key(basename(run_dirs)) %in% normalize_key(strategies)]
    }

    for (run_dir in run_dirs) {
      run_name <- basename(run_dir)
      tool_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
      tool_dirs <- tool_dirs[file.info(tool_dirs)$isdir]
      if (!is.null(tools)) {
        tool_dirs <- tool_dirs[normalize_key(basename(tool_dirs)) %in% normalize_key(tools)]
      }

      for (tool_dir in tool_dirs) {
        tool_name <- basename(tool_dir)
        res <- process_one_tool(
          sample_name = sample_name,
          run_name = run_name,
          method_specs = tool_method_specs(tool_name, input_basenames),
          tool_dir = tool_dir,
          selected_spots = selected_spots,
          spot_source = spot_source,
          overwrite = overwrite,
          dry_run = dry_run
        )
        summary_rows[[row_i]] <- res
        row_i <- row_i + 1L
      }
    }
  }

  if (length(summary_rows) == 0L) {
    log_msg("[Done] %s | No sample/run/tool combinations found.", dataset_name)
    return(NULL)
  }

  summary <- rbindlist(summary_rows, fill = TRUE)
  summary[, dataset := dataset_name]
  setcolorder(summary, c("dataset", setdiff(names(summary), "dataset")))

  log_msg(
    "[Done] %s | rows=%d | saved=%d | dry_run=%d | exists_skip=%d | input_missing=%d | parse_failed=%d | no_match=%d",
    dataset_name,
    nrow(summary),
    sum(summary$status == "saved"),
    sum(summary$status == "dry_run"),
    sum(summary$status == "exists_skip"),
    sum(summary$status == "input_missing"),
    sum(summary$status == "matrix_parse_failed"),
    sum(summary$status == "no_matched_spots")
  )

  summary
}

main <- function() {
  root_dir <- DEFAULT_ROOT_DIR
  datasets_dir <- DEFAULT_DATASETS_DIR
  datasets <- NULL
  samples <- NULL
  strategies <- NULL
  tools <- NULL
  decon_run <- DEFAULT_DECON_RUN
  decon_keep_file <- DEFAULT_DECON_KEEP_FILE
  input_basenames <- DEFAULT_INPUT_BASENAMES
  overwrite <- FALSE
  dry_run <- FALSE
  histo_fallback_to_decon <- FALSE

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L) {
    i <- 1L
    while (i <= length(args)) {
      key <- args[i]
      if (key == "--root_dir" && i + 1L <= length(args)) {
        root_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--datasets_dir" && i + 1L <= length(args)) {
        datasets_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--datasets" && i + 1L <= length(args)) {
        datasets <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--samples" && i + 1L <= length(args)) {
        samples <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--strategies" && i + 1L <= length(args)) {
        strategies <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--tools" && i + 1L <= length(args)) {
        tools <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--decon_run" && i + 1L <= length(args)) {
        decon_run <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--decon_keep_file" && i + 1L <= length(args)) {
        decon_keep_file <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--input_basenames" && i + 1L <= length(args)) {
        input_basenames <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--overwrite") {
        overwrite <- TRUE; i <- i + 1L; next
      }
      if (key == "--dry_run") {
        dry_run <- TRUE; i <- i + 1L; next
      }
      if (key == "--histo_fallback_to_decon") {
        histo_fallback_to_decon <- TRUE; i <- i + 1L; next
      }
      i <- i + 1L
    }
  }

  if (length(input_basenames) == 0L) stop("No input basenames provided.")

  roots <- if (!is.null(root_dir) && nzchar(root_dir)) {
    data.table(dataset = dataset_name_from_root(root_dir), root_dir = normalizePath(root_dir, mustWork = TRUE))
  } else {
    discover_dataset_roots(datasets_dir = datasets_dir, datasets = datasets)
  }
  if (nrow(roots) == 0L) {
    log_msg("[Done] No dataset roots found.")
    return(invisible(NULL))
  }

  all_rows <- list()
  for (i in seq_len(nrow(roots))) {
    res <- process_root(
      dataset_name = roots$dataset[i],
      root_dir = roots$root_dir[i],
      samples = samples,
      strategies = strategies,
      tools = tools,
      decon_run = decon_run,
      decon_keep_file = decon_keep_file,
      input_basenames = input_basenames,
      overwrite = overwrite,
      dry_run = dry_run,
      histo_fallback_to_decon = histo_fallback_to_decon
    )
    if (!is.null(res)) all_rows[[length(all_rows) + 1L]] <- res
  }

  if (length(all_rows) > 0L && nrow(roots) > 1L) {
    all_dt <- rbindlist(all_rows, fill = TRUE)
    log_msg("[Done] all_datasets=%d | rows=%d | saved=%d | dry_run=%d",
            length(unique(all_dt$dataset)),
            nrow(all_dt),
            sum(all_dt$status == "saved"),
            sum(all_dt$status == "dry_run"))
  }
}

main()
