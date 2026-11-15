#!/usr/bin/env Rscript
# Shared config and helpers for extract_subclone_and_ari.R.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ------------------------------------------------------------------------------
# Fixed config
# ------------------------------------------------------------------------------
CFG <- list(
  synthetic_name = "synthetic",
  single_cell_dir = "single_cell_data",
  spatial_dir = "spatial_data",
  sc_meta_file_name = "sc-meta_data.tsv",
  cell_by_spot_file_name = "cell_by_spot.txt",
  pred_out_file_name = "subclone_spots_x_subclones.csv",
  gt_hard_out_file = "subclone_groundtruth_spot_hard.csv",
  gt_soft_out_file = "subclone_groundtruth_spot_soft_long.csv",
  pred_tumor_only = TRUE,
  pred_tumor_prop_cutoff = 0.5,
  gt_tumor_only = TRUE,
  major_subclone_min_cell_frac = 0.05,
  ari_min_overlap = 10L,
  pred_summary_file = "subclone_matrix_summary.csv",
  gt_summary_file = "subclone_groundtruth_summary.csv",
  ari_summary_file = "subclone_ari_summary.csv"
)

# ------------------------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------------------------
count_subclones <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(0L)
  uniqueN(df$subclone)
}

subset_spots <- function(df, keep_spots) {
  if (is.null(df) || is.null(keep_spots)) return(df)
  df[df$spot_id %in% keep_spots, c("spot_id", "subclone"), drop = FALSE]
}

write_subclone_table <- function(df, out_csv) {
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  out <- df[, c("spot_id", "subclone"), drop = FALSE]
  fwrite(out, out_csv)
  TRUE
}

# ------------------------------------------------------------------------------
# Monoclone fallback helpers
#
# When a tool produces no subclone split we distinguish two causes:
#   - "ran but found 1 clone" (or crashed for lack of subclonal contrast)
#       -> inject a single-cluster prediction so ARI is scored honestly
#   - "found no tumor at all" -> genuine detection failure (missing_prediction)
# ------------------------------------------------------------------------------
make_monoclone_df <- function(spots) {
  spots <- as.character(spots)
  spots <- spots[!is.na(spots) & nzchar(spots)]
  spots <- unique(spots)
  if (length(spots) == 0) return(NULL)
  data.frame(spot_id = spots, subclone = "clone1", stringsAsFactors = FALSE)
}

# SCEVAN no-output classifier.
# Returns "monoclone" if SCEVAN ran and called tumor spots but produced no
# subclone split; "no_tumor" if SCEVAN classified everything as normal;
# "missing" if SCEVAN did not run / output unreadable.
classify_scevan_nooutput <- function(scevan_dir) {
  if (!dir.exists(scevan_dir)) return("missing")
  rds_file <- file.path(scevan_dir, "scevan_results.rds")
  if (!file.exists(rds_file)) return("missing")

  obj <- tryCatch(readRDS(rds_file), error = function(e) NULL)
  if (is.null(obj) || !is.data.frame(obj) || !"class" %in% colnames(obj)) return("missing")

  cls <- tolower(as.character(obj$class))
  if (any(cls == "tumor", na.rm = TRUE)) "monoclone" else "no_tumor"
}

# Numbat no-output classifier.
# Returns "no_cnv_detected" when Numbat crashed with the zero-variance HMM
# error ("data are essentially constant"), indicating no CNV contrast —
# these reps are skipped (missing_prediction), not injected as monoclone.
# Returns "missing" for any other absence (Numbat did not run or unreadable).
classify_numbat_nooutput <- function(numbat_dir) {
  if (!dir.exists(numbat_dir)) return("missing")
  log_file <- file.path(numbat_dir, "log.txt")
  if (!file.exists(log_file)) return("missing")

  log_text <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
  if (grepl("data are essentially constant", log_text, fixed = TRUE)) return("no_cnv_detected")
  return("missing")
}

# ------------------------------------------------------------------------------
# Predicted subclone helpers (from tool outputs)
# ------------------------------------------------------------------------------
normalize_clone_key <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- gsub("^subclone[_ -]*", "", s)
  s <- gsub("^clone[_ -]*", "", s)
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) s <- "unknown"
  s
}

relabel_to_cloneN <- function(labels) {
  keys <- vapply(labels, normalize_clone_key, character(1))
  uniq <- unique(keys)
  num <- suppressWarnings(as.numeric(uniq))
  ord <- order(is.na(num), num, uniq)
  uniq <- uniq[ord]
  mapper <- setNames(paste0("clone", seq_along(uniq)), uniq)
  unname(mapper[keys])
}

build_subclone_table <- function(spots, labels) {
  spots <- as.character(spots)
  labels <- as.character(labels)
  keep <- !is.na(spots) & nzchar(spots) & !is.na(labels) & nzchar(labels)
  if (!any(keep)) return(NULL)

  dt <- data.table(
    spot_id = spots[keep],
    raw_label = labels[keep]
  )
  dt[, subclone := relabel_to_cloneN(raw_label)]
  dt[, spot_ord := match(spot_id, unique(spot_id))]

  votes <- dt[, .N, by = .(spot_id, subclone, spot_ord)]
  votes[, clone_idx := suppressWarnings(as.integer(sub("^clone", "", subclone)))]
  votes[is.na(clone_idx), clone_idx := .Machine$integer.max]
  setorder(votes, spot_ord, -N, clone_idx, subclone)

  out <- votes[, .SD[1], by = spot_id]
  setorder(out, spot_ord)
  as.data.frame(out[, .(spot_id, subclone)])
}

read_tumor_spots_metadata <- function(rep_dir, std_in_name, tumor_only, tumor_prop_cutoff) {
  if (!tumor_only) return(NULL)

  meta_file <- file.path(rep_dir, std_in_name, "Metadata_cellprop.csv")
  if (!file.exists(meta_file)) {
    alt <- c(
      file.path(rep_dir, "Standardized_Inputs", "Metadata_cellprop.csv"),
      file.path(rep_dir, "Standardized_Inputs_TME", "Metadata_cellprop.csv")
    )
    hit <- alt[file.exists(alt)]
    if (length(hit) > 0) meta_file <- hit[1]
  }
  if (!file.exists(meta_file)) return(NULL)

  df <- fread(meta_file, data.table = FALSE)
  if (!"Tumor_prop" %in% colnames(df) || ncol(df) < 1) return(NULL)

  spot_col <- colnames(df)[1]
  tumor_ids <- as.character(df[[spot_col]][df$Tumor_prop > tumor_prop_cutoff])
  tumor_ids[!is.na(tumor_ids) & nzchar(tumor_ids)]
}

make_unique_coord_names <- function(coords) {
  coords <- as.character(coords)
  if (!any(duplicated(coords))) return(coords)

  seen <- new.env(parent = emptyenv())
  out <- character(length(coords))

  for (i in seq_along(coords)) {
    nm <- coords[i]
    if (!exists(nm, envir = seen, inherits = FALSE)) {
      assign(nm, 0, envir = seen)
      out[i] <- nm
    } else {
      k <- get(nm, envir = seen, inherits = FALSE) + 1
      assign(nm, k, envir = seen)
      out[i] <- paste0(nm, "_", k)
    }
  }
  out
}

map_starch_coords_to_spots <- function(rep_dir, label_coords) {
  meta_file <- file.path(rep_dir, "synthetic", "spatial_data", "st-meta.parquet")
  if (!file.exists(meta_file)) return(rep(NA_character_, length(label_coords)))
  if (!requireNamespace("arrow", quietly = TRUE)) {
    warning("[STARCH] arrow package missing; cannot map coordinates to spot IDs.")
    return(rep(NA_character_, length(label_coords)))
  }

  meta <- as.data.frame(arrow::read_parquet(meta_file))
  if (!all(c("index", "x", "y") %in% names(meta))) {
    warning("[STARCH] st-meta.parquet must contain: index, x, y")
    return(rep(NA_character_, length(label_coords)))
  }

  coord <- paste0(as.integer(round(meta$x)), "x", as.integer(round(meta$y)))
  coord <- make_unique_coord_names(coord)

  normalize_coord <- function(v) gsub("\\.0", "", as.character(v))
  coord_norm <- normalize_coord(coord)
  label_norm <- normalize_coord(label_coords)

  meta$index[match(label_norm, coord_norm)]
}

resolve_infercnv_dir <- function(tool_root) {
  if (!dir.exists(tool_root)) return(NULL)
  dirs <- list.dirs(tool_root, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0) return(NULL)

  infer_dirs <- dirs[grepl("^InferCNV", basename(dirs))]
  if (length(infer_dirs) == 0) return(NULL)

  exact <- infer_dirs[basename(infer_dirs) == "InferCNV"]
  if (length(exact) > 0) return(exact[1])

  no_hmm <- infer_dirs[basename(infer_dirs) == "InferCNV_noHMM"]
  if (length(no_hmm) > 0) return(no_hmm[1])

  sort(infer_dirs)[1]
}

extract_infercnv_table <- function(infer_dir) {
  if (is.null(infer_dir) || !dir.exists(infer_dir)) return(NULL)

  infer_file <- file.path(infer_dir, "infercnv_subclusters.observation_groupings.txt")
  if (!file.exists(infer_file)) return(NULL)

  # InferCNV stores spot IDs as row names in this file, and "Annotation Group"
  # is the subclone assignment we want to compare against ground truth.
  df <- read.table(
    infer_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "\"",
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0 || !"Annotation Group" %in% colnames(df)) return(NULL)

  # Some batches prefix spot IDs with the sample name (e.g. "Rep_1__Spot_500");
  # strip any leading "<prefix>__" so IDs match the ground-truth "Spot_N" form.
  spot_ids <- sub("^.*__", "", rownames(df))

  build_subclone_table(spot_ids, df[["Annotation Group"]])
}

extract_calicost_table <- function(calicost_dir) {
  if (is.null(calicost_dir) || !dir.exists(calicost_dir)) return(NULL)

  init_dirs <- list.dirs(calicost_dir, full.names = TRUE, recursive = FALSE)
  init_dirs <- init_dirs[grepl("^clone\\d+_rectangle0_w1\\.0$", basename(init_dirs))]
  if (length(init_dirs) == 0) return(NULL)

  init_dir <- init_dirs[order(file.info(init_dirs)$mtime, decreasing = TRUE)][1]
  labels_file <- file.path(init_dir, "clone_labels.tsv")
  if (!file.exists(labels_file)) return(NULL)

  dt <- fread(labels_file)
  if (ncol(dt) < 2 || nrow(dt) == 0) return(NULL)
  setnames(dt, 1:2, c("spot_id", "clone_label"))

  build_subclone_table(dt$spot_id, dt$clone_label)
}

extract_scevan_table <- function(scevan_dir) {
  if (!dir.exists(scevan_dir)) return(NULL)
  rds_file <- file.path(scevan_dir, "scevan_results.rds")
  if (!file.exists(rds_file)) return(NULL)

  obj <- readRDS(rds_file)
  if (!is.data.frame(obj) || !"subclone" %in% colnames(obj)) return(NULL)
  if (is.null(rownames(obj))) return(NULL)

  keep <- !is.na(obj$subclone)
  if (!any(keep)) return(NULL)
  build_subclone_table(rownames(obj)[keep], obj$subclone[keep])
}

extract_starch_table <- function(starch_dir, rep_dir) {
  if (!dir.exists(starch_dir)) return(NULL)
  labels_file <- file.path(starch_dir, "labels_STITCH_output.csv")
  if (!file.exists(labels_file)) return(NULL)

  labels_df <- read.csv(labels_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(labels_df) < 1 || nrow(labels_df) == 0) return(NULL)

  labels <- labels_df[[1]]
  coords <- rownames(labels_df)
  spots <- map_starch_coords_to_spots(rep_dir, coords)
  keep <- !is.na(spots) & !is.na(labels)
  if (!any(keep)) return(NULL)

  build_subclone_table(spots[keep], labels[keep])
}

extract_numbat_table <- function(numbat_dir) {
  if (!dir.exists(numbat_dir)) return(NULL)

  clone_files <- list.files(
    numbat_dir,
    pattern = "^clone_post_[0-9]+\\.tsv$",
    full.names = TRUE
  )
  if (length(clone_files) == 0) {
    clone_files <- list.files(
      numbat_dir,
      pattern = "^clone_post\\.tsv$",
      full.names = TRUE
    )
  }
  if (length(clone_files) == 0) return(NULL)

  get_idx <- function(fp) suppressWarnings(as.integer(sub(".*clone_post_([0-9]+)\\.tsv$", "\\1", basename(fp))))
  idx <- get_idx(clone_files)
  if (all(is.na(idx))) {
    clone_file <- sort(clone_files)[1]
  } else {
    clone_file <- clone_files[which.max(ifelse(is.na(idx), -Inf, idx))]
  }

  df <- fread(clone_file, data.table = FALSE)
  if (!"cell" %in% colnames(df) || nrow(df) == 0) return(NULL)

  # Keep ALL of numbat's clone_opt labels (do not pre-filter by numbat's own
  # compartment_opt tumor/normal call). Every other tool keeps all clusters and
  # lets the tumour region (pred_tumor_spots, applied downstream) define scope;
  # numbat must match that. Pre-filtering to compartment_opt=="tumor" dropped
  # tumour spots numbat happened to label "normal", collapsing genuine
  # multi-clone predictions to a single clone.

  if ("clone_opt" %in% colnames(df) && any(!is.na(df$clone_opt))) {
    return(build_subclone_table(df$cell, df$clone_opt))
  }

  p_cols <- grep("^p_[0-9]+$", colnames(df), value = TRUE)
  if (length(p_cols) == 0) return(NULL)
  p_mat <- as.matrix(df[, p_cols, drop = FALSE])
  p_mat[!is.finite(p_mat)] <- 0
  labels <- sub("^p_", "", p_cols)[max.col(p_mat, ties.method = "first")]
  build_subclone_table(df$cell, labels)
}

# ------------------------------------------------------------------------------
# Ground-truth helpers (from synthetic truth)
# ------------------------------------------------------------------------------
normalize_clone_id <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.integer(x))
}

to_clone_labels_from_ids <- function(clone_ids) {
  uniq <- sort(unique(clone_ids))
  mapper <- setNames(paste0("clone", seq_along(uniq)), as.character(uniq))
  unname(mapper[as.character(clone_ids)])
}

cell_id_to_int <- function(x) {
  suppressWarnings(as.integer(sub("^Cell_", "", as.character(x))))
}

read_sc_meta <- function(sc_meta_file) {
  if (!file.exists(sc_meta_file)) return(NULL)

  dt <- fread(sc_meta_file)
  if (ncol(dt) < 2) return(NULL)

  id_col <- names(dt)[1]
  setnames(dt, id_col, "cell_id")
  dt[, cell_id := as.character(cell_id)]

  if (!"clone_id" %in% names(dt)) return(NULL)
  dt[, clone_id := normalize_clone_id(clone_id)]

  if ("cell_type" %in% names(dt)) {
    dt[, cell_type := as.character(cell_type)]
  } else {
    dt[, cell_type := NA_character_]
  }

  dt[, cell_num := cell_id_to_int(cell_id)]
  dt
}

parse_cell_by_spot <- function(cell_by_spot_file) {
  if (!file.exists(cell_by_spot_file)) return(NULL)

  lines <- readLines(cell_by_spot_file, warn = FALSE)
  if (length(lines) == 0) return(NULL)

  out_list <- vector("list", length(lines))
  out_i <- 0L

  for (ln in lines) {
    if (!grepl(":", ln, fixed = TRUE)) next

    sp <- strsplit(ln, ":", fixed = TRUE)[[1]]
    if (length(sp) < 2) next

    spot_id <- trimws(sp[1])
    rhs <- trimws(paste(sp[-1], collapse = ":"))
    if (!nzchar(spot_id) || !nzchar(rhs)) next

    ids <- trimws(strsplit(rhs, ",", fixed = TRUE)[[1]])
    ids <- ids[nzchar(ids)]
    if (length(ids) == 0) next

    if (all(grepl("^Cell_", ids))) {
      cell_ids <- ids
    } else {
      cell_ids <- paste0("Cell_", ids)
    }

    out_i <- out_i + 1L
    out_list[[out_i]] <- data.table(
      spot_id = spot_id,
      cell_id = cell_ids
    )
  }

  if (out_i == 0L) return(NULL)
  rbindlist(out_list[seq_len(out_i)], use.names = TRUE)
}

read_lineage_table <- function(lineage_file) {
  if (!file.exists(lineage_file)) return(NULL)

  dt <- fread(lineage_file)
  if (ncol(dt) < 3) return(NULL)

  pick_col <- function(nms, candidates, fallback_idx) {
    hit <- candidates[candidates %in% nms]
    if (length(hit) > 0L) return(hit[1])
    nms[fallback_idx]
  }

  child_col <- pick_col(
    names(dt),
    c("child_row", "child_id", "child", "cell_num", "cell_id"),
    2L
  )
  parent_col <- pick_col(
    names(dt),
    c("parent_row", "parent_id", "parent", "parent_num"),
    3L
  )

  out <- dt[, .(
    cell_num = cell_id_to_int(get(child_col)),
    parent_num = cell_id_to_int(get(parent_col))
  )]
  out <- out[!is.na(cell_num)]
  out <- unique(out, by = "cell_num")
  out
}

read_founder_cells <- function(event_log_file) {
  if (!file.exists(event_log_file)) return(integer(0))

  dt <- fread(event_log_file)
  if (!all(c("iter", "cell_id") %in% names(dt))) return(integer(0))

  keep <- dt$iter > 0
  if ("cell_type" %in% names(dt)) {
    keep <- keep & tolower(as.character(dt$cell_type)) == "tumor"
  }
  if ("context" %in% names(dt)) {
    keep <- keep & !(tolower(as.character(dt$context)) %in% c("init_shared", "init_private"))
  }

  out <- suppressWarnings(as.integer(dt$cell_id[keep]))
  unique(out[!is.na(out)])
}

derive_cell_subclone_ids <- function(sc_dt, lineage_dt, founder_cells, min_frac, tumor_only = TRUE) {
  dt <- copy(sc_dt)
  if (tumor_only && "cell_type" %in% names(dt) && any(!is.na(dt$cell_type))) {
    dt <- dt[tolower(cell_type) == "tumor"]
  }
  dt <- dt[!is.na(clone_id) & !is.na(cell_num), .(cell_id, cell_num, clone_id)]

  if (nrow(dt) == 0 || is.null(lineage_dt) || nrow(lineage_dt) == 0) {
    return(list(
      cell_labels = NULL,
      founder_stats = data.table(),
      major_founders = data.table(),
      n_tumor_cells = 0L,
      n_major_founders = 0L
    ))
  }

  lineage_cells <- unique(as.integer(lineage_dt$cell_num))
  dt <- dt[cell_num %in% lineage_cells]
  if (nrow(dt) == 0) {
    return(list(
      cell_labels = NULL,
      founder_stats = data.table(),
      major_founders = data.table(),
      n_tumor_cells = 0L,
      n_major_founders = 0L
    ))
  }

  # "Survived tumor cells" are the tumor cells present in sc-meta / spot mapping.
  survived_tumor_cells <- unique(as.integer(dt$cell_num))
  n_tumor_cells <- length(survived_tumor_cells)

  # lineage_all supplies full parent chain, including dead ancestors.
  parent_map <- setNames(as.integer(lineage_dt$parent_num), as.character(lineage_dt$cell_num))

  get_ancestors <- local({
    cache <- new.env(hash = TRUE, parent = emptyenv())
    function(cell_num) {
      key <- as.character(cell_num)
      if (exists(key, envir = cache, inherits = FALSE)) {
        return(get(key, envir = cache, inherits = FALSE))
      }
      out <- integer(0)
      cur <- as.integer(cell_num)
      seen <- new.env(hash = TRUE, parent = emptyenv())
      while (!is.na(cur) && cur >= 0) {
        cur_key <- as.character(cur)
        if (exists(cur_key, envir = seen, inherits = FALSE)) break
        assign(cur_key, TRUE, envir = seen)
        out <- c(out, cur)
        p <- parent_map[cur_key]
        if (length(p) == 0 || is.na(p) || p < 0) break
        cur <- as.integer(p)
      }
      assign(key, out, envir = cache)
      out
    }
  })

  founder_cells <- unique(as.integer(founder_cells))
  founder_cells <- founder_cells[!is.na(founder_cells) & founder_cells %in% lineage_cells]

  founder_stats <- data.table(founder_cell = founder_cells)
  founder_stats[, tumor_desc := 0L]
  founder_stats[, frac := 0.0]

  if (nrow(founder_stats) > 0 && n_tumor_cells > 0) {
    founder_counts <- integer(length(founder_cells))
    names(founder_counts) <- as.character(founder_cells)

    for (cid in survived_tumor_cells) {
      anc <- get_ancestors(cid)
      if (length(anc) == 0) next
      hit <- intersect(anc, founder_cells)
      if (length(hit) > 0) {
        founder_counts[as.character(hit)] <- founder_counts[as.character(hit)] + 1L
      }
    }

    founder_stats[, tumor_desc := as.integer(founder_counts[as.character(founder_cell)])]
    founder_stats[is.na(tumor_desc), tumor_desc := 0L]
    founder_stats[, frac := tumor_desc / n_tumor_cells]
  }

  major_founders <- founder_stats[frac >= min_frac]
  if (nrow(major_founders) > 0) {
    node_depth <- function(node) {
      anc <- get_ancestors(node)
      if (length(anc) == 0) return(0L)
      as.integer(length(anc) - 1L)
    }

    major_founders[, depth := vapply(founder_cell, node_depth, integer(1))]
    setorder(major_founders, depth, founder_cell)

    seed_max <- suppressWarnings(max(dt$clone_id, na.rm = TRUE))
    if (!is.finite(seed_max)) seed_max <- 0L
    major_founders[, new_clone_id := as.integer(seed_max + seq_len(.N))]
  } else {
    major_founders <- founder_stats[0]
    major_founders[, depth := integer(0)]
    major_founders[, new_clone_id := integer(0)]
  }

  founder_new_map <- setNames(as.integer(major_founders$new_clone_id), as.character(major_founders$founder_cell))
  assign_one <- function(cell_num, base_clone) {
    anc <- get_ancestors(cell_num)
    if (length(anc) == 0) return(as.integer(base_clone))
    for (cur in anc) {
      hit <- founder_new_map[as.character(cur)]
      if (length(hit) > 0 && !is.na(hit)) return(as.integer(hit))
    }
    as.integer(base_clone)
  }

  dt[, subclone_id := mapply(assign_one, cell_num, clone_id, USE.NAMES = FALSE)]
  list(
    cell_labels = dt[, .(cell_id, subclone_id)],
    founder_stats = founder_stats,
    major_founders = major_founders,
    n_tumor_cells = as.integer(n_tumor_cells),
    n_major_founders = as.integer(nrow(major_founders))
  )
}

extract_gt_tables <- function(sc_dt, spot_cell_dt, lineage_dt, founder_cells, min_frac = 0.05, tumor_only = TRUE) {
  assigns <- derive_cell_subclone_ids(
    sc_dt = sc_dt,
    lineage_dt = lineage_dt,
    founder_cells = founder_cells,
    min_frac = min_frac,
    tumor_only = tumor_only
  )

  if (is.null(assigns$cell_labels) || nrow(assigns$cell_labels) == 0) {
    return(list(
      hard = NULL,
      soft = NULL,
      n_tumor_cells = assigns$n_tumor_cells,
      n_major_founders = assigns$n_major_founders
    ))
  }

  merged <- merge(
    spot_cell_dt,
    assigns$cell_labels,
    by = "cell_id",
    all = FALSE
  )
  if (nrow(merged) == 0) {
    return(list(
      hard = NULL,
      soft = NULL,
      n_tumor_cells = assigns$n_tumor_cells,
      n_major_founders = assigns$n_major_founders
    ))
  }

  counts <- merged[, .N, by = .(spot_id, subclone_id)]
  totals <- counts[, .(tumor_cells = sum(N)), by = spot_id]
  counts <- merge(counts, totals, by = "spot_id", all.x = TRUE)
  counts[, prop := N / tumor_cells]
  counts[, subclone := to_clone_labels_from_ids(subclone_id)]

  setorder(counts, spot_id, -N, subclone_id)
  hard <- counts[, .SD[1], by = spot_id][, .(spot_id, subclone)]
  soft <- counts[, .(spot_id, subclone, n_cells = N, tumor_cells, prop)]
  setorder(soft, spot_id, subclone)

  list(
    hard = hard,
    soft = soft,
    n_tumor_cells = assigns$n_tumor_cells,
    n_major_founders = assigns$n_major_founders
  )
}

# ------------------------------------------------------------------------------
# ARI helpers
# ------------------------------------------------------------------------------
choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(labels_true, labels_pred) {
  labels_true <- as.character(labels_true)
  labels_pred <- as.character(labels_pred)
  keep <- !is.na(labels_true) & !is.na(labels_pred) & nzchar(labels_true) & nzchar(labels_pred)
  labels_true <- labels_true[keep]
  labels_pred <- labels_pred[keep]

  n <- length(labels_true)
  if (n < 2) return(NA_real_)

  tab <- table(labels_true, labels_pred)
  nij <- as.numeric(tab)
  ai <- rowSums(tab)
  bj <- colSums(tab)

  sum_nij <- sum(choose2(nij))
  sum_ai <- sum(choose2(ai))
  sum_bj <- sum(choose2(bj))
  total_pairs <- choose2(n)
  if (total_pairs <= 0) return(NA_real_)

  expected <- (sum_ai * sum_bj) / total_pairs
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  numer <- sum_nij - expected

  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
    if (abs(numer) < .Machine$double.eps) return(1.0)
    return(NA_real_)
  }
  numer / denom
}

compute_ari_row <- function(rep_name, tool, gt_df, pred_df, gt_file, pred_file, min_overlap) {
  n_gt <- if (is.null(gt_df)) 0L else nrow(gt_df)
  n_pred <- if (is.null(pred_df)) 0L else nrow(pred_df)
  n_gt_sub <- count_subclones(gt_df)
  n_pred_sub <- count_subclones(pred_df)

  base_row <- data.frame(
    rep = rep_name,
    tool = tool,
    gt_file = ifelse(is.null(gt_file), NA, gt_file),
    pred_file = ifelse(is.null(pred_file), NA, pred_file),
    n_gt_spots = n_gt,
    n_pred_spots = n_pred,
    n_gt_subclones = n_gt_sub,
    n_pred_subclones = n_pred_sub,
    n_overlap_spots = 0L,
    gt_coverage = NA_real_,
    pred_coverage = NA_real_,
    ari = NA_real_,
    status = "missing_data",
    stringsAsFactors = FALSE
  )

  if (is.null(gt_df) || n_gt == 0) {
    base_row$status <- "missing_ground_truth"
    return(base_row)
  }
  if (is.null(pred_df) || n_pred == 0) {
    base_row$status <- "missing_prediction"
    return(base_row)
  }

  m <- merge(
    gt_df[, c("spot_id", "subclone"), drop = FALSE],
    pred_df[, c("spot_id", "subclone"), drop = FALSE],
    by = "spot_id",
    suffixes = c("_gt", "_pred")
  )
  n_overlap <- nrow(m)
  base_row$n_overlap_spots <- n_overlap
  base_row$gt_coverage <- ifelse(n_gt > 0, n_overlap / n_gt, NA_real_)
  base_row$pred_coverage <- ifelse(n_pred > 0, n_overlap / n_pred, NA_real_)

  if (n_overlap < min_overlap) {
    base_row$status <- "insufficient_overlap"
    return(base_row)
  }

  ari_val <- adjusted_rand_index(m$subclone_gt, m$subclone_pred)
  base_row$ari <- ari_val
  base_row$status <- ifelse(is.na(ari_val), "ari_na", "ok")
  base_row
}
