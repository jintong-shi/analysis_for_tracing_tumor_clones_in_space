#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

DEFAULT_RESULTS_DIR <- "<DATA_ROOT>/realdata/results/task2_subclone_vs_spatial"
DEFAULT_SUMMARY_FILE <- file.path(DEFAULT_RESULTS_DIR, "subclone_predictions_summary.csv")
DEFAULT_OUT_DIR <- DEFAULT_RESULTS_DIR
DEFAULT_OUT_METRIC <- "subclone_spatial_assortativity_summary.csv"
DEFAULT_OUT_NULL <- "subclone_spatial_assortativity_null.csv"
DEFAULT_N_PERM <- 1000L
DEFAULT_K <- 6L
DEFAULT_SEED <- 123L
DEFAULT_SAVE_NULL <- TRUE
DEFAULT_RESTRICT_TO_TUMOR <- TRUE
TUMOR_MASK_PREF_ORDER <- c("HistoTumor_Benign-ref", "Decon_StromaImm-ref")

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

parse_bool_arg <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  y <- tolower(trimws(x))
  if (y %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (y %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

normalize_spot_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("\\.", "-", x)
}

log_msg <- function(...) cat(sprintf(...), "\n")

find_coord_file <- function(run_dir) {
  sample_dir <- dirname(run_dir)
  allspots_dir <- file.path(sample_dir, "AllSpots_no-ref")

  candidates <- c(
    file.path(run_dir, "STARCH", "starch_10x_subset", "tissue_positions_list.csv"),
    file.path(run_dir, "STARCH", "tissue_positions_list.csv"),
    file.path(allspots_dir, "STARCH", "starch_10x_subset", "tissue_positions_list.csv"),
    file.path(allspots_dir, "STARCH", "tissue_positions_list.csv")
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) > 0L) return(candidates[1])

  fallback <- list.files(
    sample_dir,
    pattern = "^tissue_positions_list\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(fallback) > 0L) return(sort(fallback)[1])
  NA_character_
}

read_barcodes_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  bc <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  bc <- normalize_spot_ids(bc)
  bc[nzchar(bc)]
}

find_tumor_mask <- function(sample_dir) {
  for (rn in TUMOR_MASK_PREF_ORDER) {
    rdir <- file.path(sample_dir, rn)
    keep_f <- file.path(rdir, "keep_barcodes.txt")
    ref_f  <- file.path(rdir, "ref_barcodes.txt")
    if (file.exists(keep_f) && file.exists(ref_f)) {
      keep <- read_barcodes_file(keep_f)
      ref  <- read_barcodes_file(ref_f)
      mask <- setdiff(unique(keep), unique(ref))
      if (length(mask) > 0L) return(list(spots = mask, source = rn))
    }
  }
  list(spots = character(0), source = NA_character_)
}

read_tissue_positions <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)

  dt <- tryCatch(
    fread(path, header = FALSE, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L || ncol(dt) < 4L) return(NULL)

  first_val <- tolower(trimws(as.character(dt[1, 1])))
  if (first_val %in% c("barcode", "barcodes")) {
    dt <- dt[-1, , drop = FALSE]
  }
  if (nrow(dt) == 0L) return(NULL)

  spot_id <- normalize_spot_ids(dt[[1]])
  array_row <- suppressWarnings(as.numeric(dt[[3]]))
  array_col <- suppressWarnings(as.numeric(dt[[4]]))
  pix_row <- if (ncol(dt) >= 5L) suppressWarnings(as.numeric(dt[[5]])) else rep(NA_real_, nrow(dt))
  pix_col <- if (ncol(dt) >= 6L) suppressWarnings(as.numeric(dt[[6]])) else rep(NA_real_, nrow(dt))

  # Prefer array coordinates; fallback to pixel coordinates.
  x <- ifelse(is.finite(array_col), array_col, pix_col)
  y <- ifelse(is.finite(array_row), array_row, pix_row)

  out <- data.table(
    spot_id = spot_id,
    x = x,
    y = y
  )
  out <- out[!is.na(spot_id) & nzchar(spot_id) & is.finite(x) & is.finite(y)]
  if (nrow(out) == 0L) return(NULL)
  unique(out, by = "spot_id")
}

read_subclone_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)

  dt <- tryCatch(
    fread(path, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L || ncol(dt) < 2L) return(NULL)

  nms <- colnames(dt)
  key <- tolower(gsub("[^a-z0-9]+", "", nms))
  spot_idx <- which(key %in% c("spotid", "barcode", "spot", "cell"))
  sub_idx <- which(key %in% c("subclone", "clone", "cluster", "label"))

  if (length(spot_idx) == 0L) spot_idx <- 1L
  if (length(sub_idx) == 0L) sub_idx <- setdiff(seq_len(ncol(dt)), spot_idx[1])[1]
  if (is.na(sub_idx) || length(sub_idx) == 0L) return(NULL)

  out <- data.table(
    spot_id = normalize_spot_ids(dt[[spot_idx[1]]]),
    subclone = as.character(dt[[sub_idx[1]]])
  )
  out <- out[!is.na(spot_id) & nzchar(spot_id) & !is.na(subclone) & nzchar(subclone)]
  if (nrow(out) == 0L) return(NULL)
  unique(out, by = "spot_id")
}

build_knn_edges <- function(coords_dt, k = 6L) {
  n <- nrow(coords_dt)
  if (is.null(coords_dt) || n < 2L) return(matrix(integer(0), ncol = 2L))

  k_eff <- as.integer(min(max(1L, k), n - 1L))
  xy <- as.matrix(coords_dt[, .(x, y)])

  nbr_idx <- NULL
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = xy, k = k_eff + 1L)
    nbr_idx <- nn$nn.idx[, -1L, drop = FALSE]
  } else {
    d <- as.matrix(stats::dist(xy))
    diag(d) <- Inf
    nbr_idx <- t(apply(d, 1, function(v) order(v)[seq_len(k_eff)]))
    if (k_eff == 1L) nbr_idx <- matrix(nbr_idx, ncol = 1L)
  }

  from <- rep(seq_len(n), each = k_eff)
  to <- as.integer(t(nbr_idx))
  edges <- cbind(pmin(from, to), pmax(from, to))
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  if (nrow(edges) == 0L) return(matrix(integer(0), ncol = 2L))

  edges <- unique(edges)
  edges[order(edges[, 1], edges[, 2]), , drop = FALSE]
}

assortativity_nominal_undirected <- function(labels_int, edges) {
  if (is.null(edges) || nrow(edges) == 0L) return(NA_real_)
  if (length(labels_int) == 0L || all(is.na(labels_int))) return(NA_real_)

  u <- labels_int[edges[, 1]]
  v <- labels_int[edges[, 2]]
  keep <- is.finite(u) & is.finite(v)
  if (!any(keep)) return(NA_real_)
  u <- as.integer(u[keep])
  v <- as.integer(v[keep])

  k <- max(labels_int, na.rm = TRUE)
  if (!is.finite(k) || k < 2L) return(NA_real_)

  # Build directed edge-mixing matrix by counting both (u,v) and (v,u).
  bins_uv <- (u - 1L) * k + v
  bins_vu <- (v - 1L) * k + u
  counts <- tabulate(bins_uv, nbins = k * k) + tabulate(bins_vu, nbins = k * k)
  e <- matrix(counts, nrow = k, ncol = k, byrow = TRUE)
  e <- e / sum(e)

  a <- rowSums(e)
  q <- sum(a * a)
  if (!is.finite(q) || abs(1 - q) < 1e-12) return(NA_real_)

  (sum(diag(e)) - q) / (1 - q)
}

permute_assortativity <- function(labels_int, edges, n_perm = 1000L) {
  obs <- assortativity_nominal_undirected(labels_int, edges)
  if (!is.finite(obs)) {
    return(list(
      obs = obs,
      null = rep(NA_real_, n_perm),
      null_mean = NA_real_,
      null_sd = NA_real_,
      z = NA_real_,
      p_two_sided = NA_real_
    ))
  }

  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null[i] <- assortativity_nominal_undirected(sample(labels_int, replace = FALSE), edges)
  }

  mu <- mean(null, na.rm = TRUE)
  sdv <- stats::sd(null, na.rm = TRUE)
  z <- if (is.finite(sdv) && sdv > 0) (obs - mu) / sdv else NA_real_

  valid_n <- sum(is.finite(null))
  p_two <- if (valid_n > 0L) {
    (sum(abs(null) >= abs(obs), na.rm = TRUE) + 1) / (valid_n + 1)
  } else {
    NA_real_
  }

  list(
    obs = obs,
    null = null,
    null_mean = mu,
    null_sd = sdv,
    z = z,
    p_two_sided = p_two
  )
}

main <- function() {
  summary_file <- DEFAULT_SUMMARY_FILE
  out_dir <- DEFAULT_OUT_DIR
  out_metric <- DEFAULT_OUT_METRIC
  out_null <- DEFAULT_OUT_NULL
  n_perm <- DEFAULT_N_PERM
  k_neighbors <- DEFAULT_K
  seed <- DEFAULT_SEED
  save_null <- DEFAULT_SAVE_NULL
  restrict_to_tumor <- DEFAULT_RESTRICT_TO_TUMOR
  tools <- NULL
  samples <- NULL
  runs <- NULL
  datasets <- NULL

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L) {
    i <- 1L
    while (i <= length(args)) {
      key <- args[i]
      if (key == "--summary_file" && i + 1L <= length(args)) {
        summary_file <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_dir" && i + 1L <= length(args)) {
        out_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_metric" && i + 1L <= length(args)) {
        out_metric <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--out_null" && i + 1L <= length(args)) {
        out_null <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--n_permutation" && i + 1L <= length(args)) {
        n_perm <- as.integer(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--k_neighbors" && i + 1L <= length(args)) {
        k_neighbors <- as.integer(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--seed" && i + 1L <= length(args)) {
        seed <- as.integer(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--save_null" && i + 1L <= length(args)) {
        save_null <- parse_bool_arg(args[i + 1L], default = DEFAULT_SAVE_NULL); i <- i + 2L; next
      }
      if (key == "--restrict_to_tumor" && i + 1L <= length(args)) {
        restrict_to_tumor <- parse_bool_arg(args[i + 1L], default = DEFAULT_RESTRICT_TO_TUMOR); i <- i + 2L; next
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
      if (key == "--datasets" && i + 1L <= length(args)) {
        datasets <- parse_csv_arg(args[i + 1L]); i <- i + 2L; next
      }
      if (key %in% c("--help", "-h")) {
        cat(
          "Usage: Rscript 07_evaluate_subclone_spatial_assortativity.R [options]\n",
          "Options:\n",
          "  --summary_file <path>  Input summary from script 06\n",
          "  --out_dir <path>       Output directory\n",
          "  --out_metric <name>    Output metric CSV filename\n",
          "  --out_null <name>      Output null CSV filename\n",
          "  --n_permutation <int>  Number of permutations [default 1000]\n",
          "  --k_neighbors <int>    k in kNN graph [default 6]\n",
          "  --seed <int>           Random seed [default 123]\n",
          "  --save_null <bool>     Save null distribution CSV [default TRUE]\n",
          "  --restrict_to_tumor <bool>  Restrict assortativity to per-sample tumor mask [default TRUE]\n",
          "  --tools <csv>          Optional tool whitelist\n",
          "  --samples <csv>        Optional sample whitelist\n",
          "  --runs <csv>           Optional run whitelist\n",
          "  --datasets <csv>       Optional dataset whitelist (when summary has dataset col)\n",
          sep = ""
        )
        quit(save = "no", status = 0L)
      }
      stop("Unknown or incomplete argument: ", key)
    }
  }

  if (!file.exists(summary_file)) stop("summary_file not found: ", summary_file)
  if (!is.finite(n_perm) || n_perm < 1L) stop("n_permutation must be >= 1")
  if (!is.finite(k_neighbors) || k_neighbors < 1L) stop("k_neighbors must be >= 1")

  set.seed(seed)

  sm <- fread(summary_file)
  if (!all(c("sample", "run", "tool", "run_dir") %in% names(sm))) {
    stop("summary_file missing required columns: sample, run, tool, run_dir")
  }
  if (!"pred_file" %in% names(sm)) {
    if (!"tool_dir" %in% names(sm)) {
      stop("summary_file must contain pred_file or tool_dir")
    }
    sm[, pred_file := file.path(tool_dir, "subclone_spots_x_subclones.csv")]
  }
  if (!"dataset" %in% names(sm)) {
    sm[, dataset := NA_character_]
  }

  sm <- sm[status == "ok"]
  if (!is.null(tools)) sm <- sm[tool %in% tools]
  if (!is.null(samples)) sm <- sm[sample %in% samples]
  if (!is.null(runs)) sm <- sm[run %in% runs]
  if (!is.null(datasets)) sm <- sm[dataset %in% datasets]

  if (nrow(sm) == 0L) {
    stop("No valid rows to evaluate after filtering.")
  }
  setorder(sm, dataset, sample, run, tool)

  metric_rows <- vector("list", nrow(sm))
  null_rows <- list()
  tumor_mask_cache <- list()

  for (i in seq_len(nrow(sm))) {
    row <- sm[i]
    dataset_id <- as.character(row$dataset)
    sample_id <- as.character(row$sample)
    run_id <- as.character(row$run)
    tool_id <- as.character(row$tool)
    run_dir <- as.character(row$run_dir)
    pred_file <- as.character(row$pred_file)

    log_msg("[%03d/%03d] %s / %s / %s / %s",
            i, nrow(sm), dataset_id, sample_id, run_id, tool_id)

    coord_file <- find_coord_file(run_dir)
    pred <- read_subclone_file(pred_file)
    coord <- read_tissue_positions(coord_file)

    tumor_source <- NA_character_
    n_tumor_mask <- 0L
    tumor_spots <- character(0)
    if (restrict_to_tumor) {
      sample_dir <- dirname(run_dir)
      cache_key <- paste(dataset_id, sample_id, sep = "::")
      if (is.null(tumor_mask_cache[[cache_key]])) {
        tumor_mask_cache[[cache_key]] <- find_tumor_mask(sample_dir)
      }
      tumor_info <- tumor_mask_cache[[cache_key]]
      tumor_spots <- tumor_info$spots
      tumor_source <- tumor_info$source
      n_tumor_mask <- length(tumor_spots)
    }

    status <- "ok"
    detail <- "assortativity_computed"
    n_pred <- if (is.null(pred)) 0L else nrow(pred)
    n_coord <- if (is.null(coord)) 0L else nrow(coord)
    n_used <- 0L
    n_edges <- 0L
    n_sub <- 0L
    obs <- NA_real_
    mu <- NA_real_
    sdv <- NA_real_
    z <- NA_real_
    p_two <- NA_real_

    if (is.null(pred) || nrow(pred) == 0L) {
      status <- "pred_missing_or_empty"
      detail <- "prediction_file_empty_or_invalid"
    } else if (is.null(coord) || nrow(coord) == 0L) {
      status <- "coord_missing_or_empty"
      detail <- "tissue_positions_missing_or_invalid"
    } else if (restrict_to_tumor && n_tumor_mask == 0L) {
      status <- "tumor_mask_missing"
      detail <- "no_histotumor_or_decon_stromaimm_for_sample"
    } else {
      dt <- merge(pred, coord, by = "spot_id")
      if (restrict_to_tumor) {
        dt <- dt[spot_id %in% tumor_spots]
      }
      n_used <- nrow(dt)
      n_sub <- uniqueN(dt$subclone)

      if (n_used < 10L) {
        status <- "too_few_spots"
        detail <- "need_at_least_10_spots_after_merge"
      } else if (n_sub < 2L) {
        status <- "single_subclone"
        detail <- "need_at_least_2_subclones_for_assortativity"
      } else {
        setorder(dt, spot_id)
        edges <- build_knn_edges(dt[, .(x, y)], k = k_neighbors)
        n_edges <- nrow(edges)
        if (n_edges < 1L) {
          status <- "graph_empty"
          detail <- "no_edges_in_knn_graph"
        } else {
          labels_int <- as.integer(factor(dt$subclone, levels = sort(unique(dt$subclone))))
          out <- permute_assortativity(labels_int, edges, n_perm = n_perm)
          obs <- out$obs
          mu <- out$null_mean
          sdv <- out$null_sd
          z <- out$z
          p_two <- out$p_two_sided

          if (save_null) {
            null_rows[[length(null_rows) + 1L]] <- data.table(
              dataset = dataset_id,
              sample = sample_id,
              run = run_id,
              tool = tool_id,
              permutation_id = seq_len(n_perm),
              assortativity_null = out$null
            )
          }
        }
      }
    }

    metric_rows[[i]] <- data.table(
      dataset = dataset_id,
      sample = sample_id,
      run = run_id,
      tool = tool_id,
      run_dir = run_dir,
      pred_file = pred_file,
      coord_file = coord_file,
      restrict_to_tumor = isTRUE(restrict_to_tumor),
      tumor_mask_source = tumor_source,
      n_tumor_mask = as.integer(n_tumor_mask),
      n_spots_pred = as.integer(n_pred),
      n_spots_coord = as.integer(n_coord),
      n_spots_used = as.integer(n_used),
      n_subclones = as.integer(n_sub),
      n_edges = as.integer(n_edges),
      k_neighbors = as.integer(k_neighbors),
      n_permutation = as.integer(n_perm),
      assortativity_obs = obs,
      null_mean = mu,
      null_sd = sdv,
      z_score = z,
      p_value_two_sided = p_two,
      status = status,
      detail = detail
    )
  }

  metric_dt <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
  setorder(metric_dt, dataset, sample, run, tool)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  metric_path <- file.path(out_dir, out_metric)
  fwrite(metric_dt, metric_path)
  log_msg("[Done] Metric table: %s", metric_path)

  if (save_null) {
    null_dt <- if (length(null_rows) > 0L) {
      rbindlist(null_rows, use.names = TRUE, fill = TRUE)
    } else {
      data.table(
        dataset = character(),
        sample = character(),
        run = character(),
        tool = character(),
        permutation_id = integer(),
        assortativity_null = numeric()
      )
    }
    setorder(null_dt, dataset, sample, run, tool, permutation_id)
    null_path <- file.path(out_dir, out_null)
    fwrite(null_dt, null_path)
    log_msg("[Done] Null table: %s", null_path)
  }
}

main()

