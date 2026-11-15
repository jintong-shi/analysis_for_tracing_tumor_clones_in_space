#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

normalize_spot_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("\\.", "-", x)
}

count_subclones <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(0L)
  uniqueN(df$subclone)
}

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
  spots <- normalize_spot_ids(spots)
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

resolve_tool_dir_by_aliases <- function(run_dir, aliases) {
  if (!dir.exists(run_dir)) return(NA_character_)
  dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0L) return(NA_character_)

  dir_key <- normalize_key(basename(dirs))
  alias_key <- normalize_key(aliases)

  hit <- dirs[dir_key %in% alias_key]
  if (length(hit) > 0L) return(sort(hit)[1])

  starts <- dirs[vapply(dir_key, function(k) any(startsWith(k, alias_key)), logical(1))]
  if (length(starts) > 0L) return(sort(starts)[1])

  NA_character_
}

resolve_infercnv_dir <- function(run_dir) {
  if (!dir.exists(run_dir)) return(NA_character_)
  dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0L) return(NA_character_)

  base <- basename(dirs)
  key <- normalize_key(base)
  infer_dirs <- dirs[grepl("^infercnv", key)]
  if (length(infer_dirs) == 0L) return(NA_character_)

  infer_base <- basename(infer_dirs)
  infer_key <- normalize_key(infer_base)
  exact <- infer_dirs[infer_key == "infercnv"]
  if (length(exact) > 0L) return(exact[1])

  no_hmm <- infer_dirs[infer_key == "infercnvnohmm"]
  if (length(no_hmm) > 0L) return(no_hmm[1])

  sort(infer_dirs)[1]
}

resolve_infercnv_named_dir <- function(run_dir, target_key) {
  if (!dir.exists(run_dir)) return(NA_character_)
  dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0L) return(NA_character_)

  key <- normalize_key(basename(dirs))
  hit <- dirs[key == normalize_key(target_key)]
  if (length(hit) > 0L) return(hit[1])

  NA_character_
}

infercnv_grouping_files <- function(infer_dir, subcluster_only = FALSE) {
  files <- if (isTRUE(subcluster_only)) {
    "infercnv_subclusters.observation_groupings.txt"
  } else {
    c(
      "infercnv_subclusters.observation_groupings.txt",
      "infercnv.observation_groupings.txt",
      "infercnv.preliminary.observation_groupings.txt"
    )
  }
  file.path(infer_dir, files)
}

coord_key_from_label <- function(x) {
  parts <- strsplit(as.character(x), "x", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2L) return(NA_character_)
    a <- suppressWarnings(as.integer(round(as.numeric(p[1]))))
    b <- suppressWarnings(as.integer(round(as.numeric(p[2]))))
    if (is.na(a) || is.na(b)) return(NA_character_)
    paste0(a, "x", b)
  }, character(1))
}

extract_infercnv_table <- function(infer_dir, subcluster_only = FALSE) {
  if (is.na(infer_dir) || !dir.exists(infer_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  cand <- infercnv_grouping_files(infer_dir, subcluster_only = subcluster_only)
  cand <- cand[file.exists(cand)]
  if (length(cand) == 0L) {
    detail <- if (isTRUE(subcluster_only)) "subcluster_groupings_file_missing" else "groupings_file_missing"
    return(list(df = NULL, input_file = NA_character_, detail = detail))
  }

  for (f in cand) {
    raw <- tryCatch(
      read.table(
        f,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        quote = "\"",
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(raw) || nrow(raw) == 0L) next

    cands <- c("Annotation Group", "Subcluster", "cluster", "Clone", "Dendrogram Group")
    nms <- colnames(raw)
    idx <- match(normalize_key(cands), normalize_key(nms))
    idx <- idx[!is.na(idx)]
    label_col <- if (length(idx) > 0L) nms[idx[1]] else nms[ncol(raw)]

    out <- build_subclone_table(rownames(raw), raw[[label_col]])
    if (!is.null(out) && nrow(out) > 0L) {
      return(list(df = out, input_file = f, detail = paste0("parsed:", label_col)))
    }
  }

  list(df = NULL, input_file = cand[1], detail = "no_usable_rows")
}

extract_infercnv_with_regular_fallback <- function(run_dir) {
  primary_dir <- resolve_infercnv_named_dir(run_dir, "infercnv")
  regular_dir <- resolve_infercnv_named_dir(run_dir, "infercnv_regular")

  if (!is.na(primary_dir)) {
    primary_sub <- extract_infercnv_table(primary_dir, subcluster_only = TRUE)
    if (!is.null(primary_sub$df) && nrow(primary_sub$df) > 0L) {
      primary_sub$tool_dir <- primary_dir
      return(primary_sub)
    }

    if (identical(primary_sub$detail, "subcluster_groupings_file_missing") && !is.na(regular_dir)) {
      fallback <- extract_infercnv_table(regular_dir, subcluster_only = TRUE)
      if (!is.null(fallback$df) && nrow(fallback$df) > 0L) {
        fallback$detail <- paste0("fallback_infercnv_regular:", fallback$detail)
        fallback$tool_dir <- primary_dir
        return(fallback)
      }
    }

    primary_all <- extract_infercnv_table(primary_dir)
    primary_all$tool_dir <- primary_dir
    return(primary_all)
  }

  tool_dir <- resolve_infercnv_dir(run_dir)
  out <- extract_infercnv_table(tool_dir)
  out$tool_dir <- tool_dir
  out
}

extract_scevan_table <- function(scevan_dir) {
  if (is.na(scevan_dir) || !dir.exists(scevan_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  rds_files <- file.path(
    scevan_dir,
    c(
      "SCEVAN_results.rds",
      "scevan_results.rds",
      file.path("output", "SCEVAN_results.rds"),
      file.path("output", "scevan_results.rds")
    )
  )
  rds_files <- rds_files[file.exists(rds_files)]
  if (length(rds_files) == 0L) {
    return(list(df = NULL, input_file = NA_character_, detail = "results_rds_missing"))
  }

  for (f in rds_files) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj) || !is.data.frame(obj) || nrow(obj) == 0L || is.null(rownames(obj))) next
    if (!"subclone" %in% colnames(obj)) next

    keep <- !is.na(obj$subclone) & nzchar(as.character(obj$subclone))
    if (!any(keep)) next

    out <- build_subclone_table(rownames(obj)[keep], obj$subclone[keep])
    if (!is.null(out) && nrow(out) > 0L) {
      return(list(df = out, input_file = f, detail = "parsed:subclone"))
    }
  }

  list(df = NULL, input_file = rds_files[1], detail = "subclone_missing_or_empty")
}

extract_starch_table <- function(starch_dir) {
  if (is.na(starch_dir) || !dir.exists(starch_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  label_files <- file.path(
    starch_dir,
    c("labels_STARCH.csv", "labels_STARCH_output.csv", "labels_STITCH_output.csv")
  )
  label_files <- label_files[file.exists(label_files)]
  if (length(label_files) == 0L) {
    return(list(df = NULL, input_file = NA_character_, detail = "labels_file_missing"))
  }
  label_file <- label_files[1]

  pos_files <- file.path(
    starch_dir,
    c(file.path("starch_10x_subset", "tissue_positions_list.csv"), "tissue_positions_list.csv")
  )
  pos_files <- pos_files[file.exists(pos_files)]
  if (length(pos_files) == 0L) {
    return(list(df = NULL, input_file = label_file, detail = "tissue_positions_missing"))
  }
  pos_file <- pos_files[1]

  labels_df <- tryCatch(read.csv(label_file, row.names = 1, check.names = FALSE), error = function(e) NULL)
  pos_df <- tryCatch(fread(pos_file, header = FALSE, data.table = FALSE), error = function(e) NULL)
  if (is.null(labels_df) || nrow(labels_df) == 0L || ncol(labels_df) < 1L) {
    return(list(df = NULL, input_file = label_file, detail = "labels_parse_failed"))
  }
  if (is.null(pos_df) || nrow(pos_df) == 0L || ncol(pos_df) < 4L) {
    return(list(df = NULL, input_file = paste(label_file, pos_file, sep = ";"), detail = "positions_parse_failed"))
  }

  labels <- labels_df[[1]]
  coords <- rownames(labels_df)
  coord_key <- coord_key_from_label(coords)

  barcodes <- normalize_spot_ids(pos_df[[1]])
  array_row <- suppressWarnings(as.integer(round(as.numeric(pos_df[[3]]))))
  array_col <- suppressWarnings(as.integer(round(as.numeric(pos_df[[4]]))))
  pos_key <- paste0(array_row, "x", array_col)

  spots <- barcodes[match(coord_key, pos_key)]
  keep <- !is.na(spots) & nzchar(spots) & !is.na(labels) & nzchar(as.character(labels))
  if (!any(keep)) {
    return(list(df = NULL, input_file = paste(label_file, pos_file, sep = ";"), detail = "no_valid_spot_mapping"))
  }

  out <- build_subclone_table(spots[keep], labels[keep])
  if (is.null(out) || nrow(out) == 0L) {
    return(list(df = NULL, input_file = paste(label_file, pos_file, sep = ";"), detail = "no_usable_rows"))
  }

  list(df = out, input_file = paste(label_file, pos_file, sep = ";"), detail = "parsed:labels")
}

extract_numbat_table <- function(numbat_dir) {
  if (is.na(numbat_dir) || !dir.exists(numbat_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  clone_files <- list.files(numbat_dir, pattern = "^clone_post_[0-9]+\\.tsv$", full.names = TRUE)
  if (length(clone_files) == 0L) {
    clone_files <- list.files(numbat_dir, pattern = "^clone_post\\.tsv$", full.names = TRUE)
  }
  if (length(clone_files) == 0L) {
    clone_files <- list.files(
      numbat_dir,
      pattern = "^clone_post(_[0-9]+)?\\.tsv$",
      full.names = TRUE,
      recursive = TRUE
    )
  }
  if (length(clone_files) == 0L) {
    return(list(df = NULL, input_file = NA_character_, detail = "clone_post_missing"))
  }

  get_idx <- function(fp) suppressWarnings(as.integer(sub(".*clone_post_([0-9]+)\\.tsv$", "\\1", basename(fp))))
  idx <- get_idx(clone_files)
  clone_file <- if (all(is.na(idx))) {
    sort(clone_files)[1]
  } else {
    clone_files[which.max(ifelse(is.na(idx), -Inf, idx))]
  }

  df <- tryCatch(fread(clone_file, data.table = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) {
    return(list(df = NULL, input_file = clone_file, detail = "clone_post_parse_failed"))
  }

  cell_col <- c("cell", "spot", "spot_id", "barcode")
  cell_col <- cell_col[cell_col %in% colnames(df)]
  if (length(cell_col) == 0L) {
    return(list(df = NULL, input_file = clone_file, detail = "spot_column_missing"))
  }
  cell_col <- cell_col[1]

  # Keep ALL of numbat's clone_opt labels (do not pre-filter by numbat's own
  # compartment_opt tumor/normal call). Every other tool keeps all clusters and
  # lets step 07's tumour mask define the tumour region; numbat must match that.
  # Pre-filtering to compartment_opt=="tumor" dropped histologically-tumour spots
  # that numbat happened to label "normal", collapsing genuine multi-clone masks
  # to a single clone (e.g. H1_2: normal-clone spots are inside the tumour mask).

  out <- NULL
  if ("clone_opt" %in% colnames(df) && any(!is.na(df$clone_opt) & nzchar(as.character(df$clone_opt)))) {
    out <- build_subclone_table(df[[cell_col]], df$clone_opt)
  } else {
    p_cols <- grep("^p_[0-9]+$", colnames(df), value = TRUE)
    if (length(p_cols) > 0L) {
      p_mat <- as.matrix(df[, p_cols, drop = FALSE])
      p_mat[!is.finite(p_mat)] <- 0
      labels <- sub("^p_", "", p_cols)[max.col(p_mat, ties.method = "first")]
      out <- build_subclone_table(df[[cell_col]], labels)
    }
  }

  if (is.null(out) || nrow(out) == 0L) {
    return(list(df = NULL, input_file = clone_file, detail = "no_usable_rows"))
  }

  list(df = out, input_file = clone_file, detail = "parsed:clone_post")
}

safe_read_rds_auto <- function(path) {
  if (!file.exists(path)) {
    return(list(obj = NULL, error = "file_missing"))
  }

  err <- NULL
  out <- tryCatch(
    {
      if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        readRDS(gzfile(path), refhook = function(x) x)
      } else {
        readRDS(path, refhook = function(x) x)
      }
    },
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(err)) return(list(obj = NULL, error = err))
  list(obj = out, error = NA_character_)
}

extract_tuscan_table <- function(tuscan_dir) {
  if (is.na(tuscan_dir) || !dir.exists(tuscan_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  obj_files <- file.path(
    tuscan_dir,
    c(
      "tumor_obj.rds.gz",
      "tumor_obj.rds",
      "obj.rds.gz",
      "obj.rds",
      "cnv_obj.rds.gz",
      "cnv_obj.rds"
    )
  )
  obj_files <- obj_files[file.exists(obj_files)]
  if (length(obj_files) == 0L) {
    return(list(df = NULL, input_file = NA_character_, detail = "tuscan_obj_missing"))
  }

  for (f in obj_files) {
    rd <- safe_read_rds_auto(f)
    obj <- rd$obj
    if (is.null(obj)) next

    slots <- attributes(obj)
    if (is.null(slots) && typeof(obj) == "list") {
      slots <- obj
    }
    if (is.null(slots)) next

    grouped <- slots[["observation_grouped_spot_indices"]]
    cnv <- slots[["cnv.data"]]
    if (is.null(cnv)) cnv <- slots[["counts.data"]]
    cnv_spots <- normalize_spot_ids(colnames(cnv))

    if (!is.list(grouped) || length(grouped) == 0L || length(cnv_spots) == 0L) next

    n_spots <- length(cnv_spots)
    lab_names <- names(grouped)
    if (is.null(lab_names) || any(!nzchar(lab_names))) {
      lab_names <- as.character(seq_along(grouped))
    }

    parts <- vector("list", length(grouped))
    pi <- 0L
    for (j in seq_along(grouped)) {
      idx <- suppressWarnings(as.integer(grouped[[j]]))
      idx <- idx[!is.na(idx) & idx >= 1L & idx <= n_spots]
      if (length(idx) == 0L) next

      pi <- pi + 1L
      parts[[pi]] <- data.table(
        spot_id = cnv_spots[idx],
        raw_label = rep(as.character(lab_names[j]), length(idx))
      )
    }

    if (pi > 0L) {
      dt <- rbindlist(parts[seq_len(pi)], use.names = TRUE, fill = TRUE)
      out <- build_subclone_table(dt$spot_id, dt$raw_label)
      if (!is.null(out) && nrow(out) > 0L) {
        return(list(df = out, input_file = f, detail = "parsed:observation_grouped_spot_indices"))
      }
    }
  }

  # Fallback parser: annotations.txt often stores <spot_id>\t<cluster_label>.
  anno_file <- file.path(tuscan_dir, "annotations.txt")
  if (file.exists(anno_file)) {
    anno <- tryCatch(
      fread(anno_file, sep = "\t", header = FALSE, data.table = FALSE),
      error = function(e) NULL
    )
    if (!is.null(anno) && nrow(anno) > 0L && ncol(anno) >= 2L) {
      out <- build_subclone_table(anno[[1]], anno[[2]])
      if (!is.null(out) && nrow(out) > 0L) {
        return(list(df = out, input_file = anno_file, detail = "parsed:annotations_fallback"))
      }
    }
  }

  list(df = NULL, input_file = obj_files[1], detail = "no_usable_rows")
}

extract_calicost_table <- function(calicost_dir) {
  if (is.na(calicost_dir) || !dir.exists(calicost_dir)) {
    return(list(df = NULL, input_file = NA_character_, detail = "tool_dir_missing"))
  }

  label_file <- file.path(calicost_dir, "clone_labels.tsv")
  if (!file.exists(label_file)) {
    return(list(df = NULL, input_file = NA_character_, detail = "clone_labels_missing"))
  }

  dt <- tryCatch(fread(label_file, data.table = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L || ncol(dt) < 2L) {
    return(list(df = NULL, input_file = label_file, detail = "clone_labels_parse_failed"))
  }

  out <- build_subclone_table(dt[[1]], dt[[2]])
  if (is.null(out) || nrow(out) == 0L) {
    return(list(df = NULL, input_file = label_file, detail = "no_usable_rows"))
  }

  list(df = out, input_file = label_file, detail = "parsed:clone_labels")
}
