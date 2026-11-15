normalize_sample_id <- function(x) {
  trimws(as.character(x))
}

normalize_dataset_id <- function(x) {
  trimws(as.character(x))
}

find_upward <- function(start_dir, relative_path) {
  if (is.null(start_dir) || !nzchar(start_dir)) return(NULL)
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  repeat {
    candidate <- file.path(cur, relative_path)
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  NULL
}

resolve_sample_overview <- function(sample_overview = NULL,
                                    dataset_root = NULL,
                                    start_dir = NULL,
                                    required = TRUE) {
  candidates <- character()
  if (!is.null(sample_overview) && nzchar(sample_overview)) {
    candidates <- c(candidates, sample_overview)
  }
  if (!is.null(dataset_root) && nzchar(dataset_root)) {
    candidates <- c(
      candidates,
      file.path(dataset_root, "..", "..", "..", "docs", "sample_overview.csv")
    )
  }
  if (!is.null(start_dir) && nzchar(start_dir)) {
    found <- find_upward(start_dir, file.path("docs", "sample_overview.csv"))
    if (!is.null(found)) candidates <- c(candidates, found)
  }

  for (p in candidates) {
    p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (file.exists(p)) return(p)
  }

  if (required) {
    stop("Cannot find docs/sample_overview.csv. Pass --sample_overview explicitly.")
  }
  NULL
}

resolve_realdata_root <- function(realdata_root = NULL,
                                  sample_overview = NULL,
                                  start_dir = NULL,
                                  required = TRUE) {
  candidates <- character()
  if (!is.null(realdata_root) && nzchar(realdata_root)) {
    candidates <- c(candidates, realdata_root)
  }
  if (!is.null(sample_overview) && nzchar(sample_overview)) {
    so <- normalizePath(sample_overview, winslash = "/", mustWork = FALSE)
    candidates <- c(
      candidates,
      file.path(dirname(dirname(so)), "realdata_unified_by_dataset")
    )
  }
  if (!is.null(start_dir) && nzchar(start_dir)) {
    found <- find_upward(start_dir, file.path("realdata_unified_by_dataset", "datasets"))
    if (!is.null(found)) candidates <- c(candidates, dirname(found))
  }

  for (p in candidates) {
    p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(p, "datasets"))) return(p)
  }

  if (required) {
    stop("Cannot find realdata_unified_by_dataset. Pass --realdata_root explicitly.")
  }
  NULL
}

task3_results_dir <- function(realdata_root) {
  file.path(realdata_root, "results", "task3_tumor_vs_histology")
}

task3_dataset_root <- function(realdata_root, dataset_id) {
  file.path(realdata_root, "datasets", dataset_id)
}

task3_truth_csv <- function(dataset_root) {
  file.path(dataset_root, "evaluation", "task3_tumor_vs_histology",
            "truth_standardized", "tumor_normal_truth.csv")
}

normalize_task3_spot_id <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\.", "-", x)
}

choose_task3_col <- function(nms, candidates, required = TRUE, label = "") {
  nms_l <- tolower(nms)
  cand_l <- tolower(candidates)
  idx <- match(cand_l, nms_l)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    if (required) {
      stop("Cannot find ", label, " column. Tried: ",
           paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  nms[idx[1]]
}

map_task3_truth_label <- function(label) {
  s <- trimws(as.character(label))
  s_l <- tolower(s)
  is_tumor <- grepl("tumou?r|cancer|malignan|gleason|^gg[0-9]", s_l)
  is_normal <- grepl("normal|benign", s_l)
  truth_label <- ifelse(
    is_tumor & !is_normal, "Tumor",
    ifelse(is_normal & !is_tumor, "Normal", "Others")
  )
  truth_bin <- ifelse(truth_label == "Tumor", 1L,
                      ifelse(truth_label == "Normal", 0L, NA_integer_))
  data.table::data.table(truth_label = truth_label, truth_bin = truth_bin)
}

read_task3_truth_from_histology_files <- function(sample_rows, dataset_id) {
  if (is.null(sample_rows) || !nrow(sample_rows)) {
    return(data.table::data.table())
  }

  out <- list()
  for (idx in seq_len(nrow(sample_rows))) {
    sample_id <- normalize_sample_id(sample_rows$sample[idx])
    histology_file <- trimws(as.character(sample_rows$histology_file[idx]))
    if (!nzchar(histology_file) || !file.exists(histology_file)) {
      next
    }

    dt <- tryCatch(data.table::fread(histology_file), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0L) next

    c_spot <- choose_task3_col(
      names(dt),
      c("Barcode", "barcode", "spot_id", "spot", "cell", "cell_id"),
      required = TRUE,
      label = "spot"
    )
    c_hist <- choose_task3_col(
      names(dt),
      c("Final_Annotations", "final_annotations", "histology", "Histology",
        "truth_label", "label", "annotation"),
      required = FALSE,
      label = "histology"
    )
    if (is.na(c_hist)) {
      c_hist <- setdiff(names(dt), c_spot)[1]
    }
    if (is.na(c_hist) || !nzchar(c_hist)) next

    truth <- map_task3_truth_label(dt[[c_hist]])
    one <- data.table::data.table(
      dataset_id = dataset_id,
      sample = sample_id,
      spot_id = normalize_task3_spot_id(dt[[c_spot]]),
      histology_label = as.character(dt[[c_hist]]),
      truth_label = truth$truth_label,
      truth_bin = truth$truth_bin,
      source_file = normalizePath(histology_file, winslash = "/", mustWork = FALSE)
    )
    out[[length(out) + 1L]] <- one[!is.na(spot_id) & nzchar(spot_id)]
  }

  if (length(out)) {
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  } else {
    data.table::data.table(
      dataset_id = character(),
      sample = character(),
      spot_id = character(),
      histology_label = character(),
      truth_label = character(),
      truth_bin = integer(),
      source_file = character()
    )
  }
}

read_task3_truth_for_dataset <- function(truth_csv,
                                         task_samples,
                                         dataset_id) {
  dataset_id_value <- dataset_id
  selected <- data.table::copy(
    task_samples[task_samples[["dataset_id"]] == dataset_id_value]
  )
  selected[, sample := normalize_sample_id(sample)]

  truth_existing <- NULL
  if (!is.null(truth_csv) && file.exists(truth_csv)) {
    truth_existing <- data.table::fread(truth_csv)
    truth_existing[, sample := normalize_sample_id(sample)]
    truth_existing[, spot_id := normalize_task3_spot_id(spot_id)]
    truth_existing <- truth_existing[sample %in% selected$sample]
  }

  existing_samples <- if (!is.null(truth_existing) && nrow(truth_existing)) {
    unique(truth_existing$sample)
  } else {
    character()
  }
  missing_samples <- setdiff(selected$sample, existing_samples)

  truth_from_overview <- if (length(missing_samples)) {
    read_task3_truth_from_histology_files(
      selected[sample %in% missing_samples],
      dataset_id = dataset_id
    )
  } else {
    data.table::data.table()
  }

  out <- data.table::rbindlist(
    Filter(Negate(is.null), list(truth_existing, truth_from_overview)),
    use.names = TRUE,
    fill = TRUE
  )
  if (!nrow(out)) {
    stop("No truth rows available for selected samples in ", dataset_id,
         ". Checked standardized truth CSV and sample_overview histology_file paths.")
  }
  out
}

task3_run_root <- function(dataset_root, required = TRUE) {
  candidates <- file.path(dataset_root, c("run_links", "run_CNV"))
  for (p in candidates) {
    if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  if (required) {
    stop("Missing run_links/run_CNV under dataset root: ", dataset_root)
  }
  NULL
}

read_histology_sample_overview <- function(sample_overview,
                                           dataset_id = NULL,
                                           require_files = TRUE) {
  if (is.null(sample_overview) || !file.exists(sample_overview)) {
    stop("Missing sample overview: ", sample_overview)
  }

  dt <- data.table::fread(sample_overview)
  data.table::setnames(dt, sub("^\ufeff", "", names(dt)))

  req <- c("sample", "datasets", "histology_available", "histology_file")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("sample_overview.csv is missing required columns: ",
         paste(miss, collapse = ", "))
  }

  dt[, sample := normalize_sample_id(sample)]
  dt[, datasets := normalize_dataset_id(datasets)]
  dt[, histology_available := tolower(trimws(as.character(histology_available)))]
  dt[, histology_file := trimws(as.character(histology_file))]

  dt <- dt[histology_available %in% c("yes", "y", "true", "1")]
  if (!is.null(dataset_id) && nzchar(dataset_id)) {
    dt <- dt[datasets == normalize_dataset_id(dataset_id)]
  }

  if (require_files && nrow(dt)) {
    missing_files <- dt[!file.exists(histology_file), unique(histology_file)]
    if (length(missing_files)) {
      stop(
        "Histology files listed in sample_overview.csv are missing:\n",
        paste(missing_files, collapse = "\n")
      )
    }
  }

  unique(dt[, .(dataset_id = datasets, sample, histology_file)])
}

read_task3_sample_table <- function(sample_overview,
                                    realdata_root = NULL,
                                    dataset_ids = NULL,
                                    samples = NULL,
                                    require_files = TRUE) {
  dt <- read_histology_sample_overview(
    sample_overview = sample_overview,
    dataset_id = NULL,
    require_files = require_files
  )
  if (!is.null(dataset_ids) && length(dataset_ids)) {
    dataset_ids <- normalize_dataset_id(dataset_ids)
    dt <- dt[dataset_id %in% dataset_ids]
  }
  if (!is.null(samples) && length(samples)) {
    samples <- normalize_sample_id(samples)
    dt <- dt[sample %in% samples]
  }
  if (!is.null(realdata_root) && nzchar(realdata_root)) {
    dt[, dataset_root := task3_dataset_root(realdata_root, dataset_id)]
  }
  unique(dt[order(dataset_id, sample)])
}

histology_sample_ids <- function(sample_overview,
                                 dataset_id = NULL,
                                 require_files = TRUE) {
  read_histology_sample_overview(
    sample_overview = sample_overview,
    dataset_id = dataset_id,
    require_files = require_files
  )[, sample]
}

filter_dt_to_histology_samples <- function(dt,
                                           sample_overview,
                                           dataset_id = NULL,
                                           sample_col = "sample",
                                           dataset_col = "dataset_id",
                                           require_files = TRUE) {
  if (!sample_col %in% names(dt)) {
    stop("Cannot filter by histology samples; missing column: ", sample_col)
  }
  keep <- read_histology_sample_overview(
    sample_overview = sample_overview,
    dataset_id = dataset_id,
    require_files = require_files
  )
  dt[, (sample_col) := normalize_sample_id(get(sample_col))]

  if (dataset_col %in% names(dt)) {
    dt[, (dataset_col) := normalize_dataset_id(get(dataset_col))]
    merge(
      dt,
      unique(keep[, .(dataset_id, sample)]),
      by.x = c(dataset_col, sample_col),
      by.y = c("dataset_id", "sample"),
      all = FALSE,
      sort = FALSE
    )
  } else {
    dt[get(sample_col) %in% keep[, sample]]
  }
}
