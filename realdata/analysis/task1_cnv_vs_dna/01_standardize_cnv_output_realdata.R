#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

EPS <- 1e-12
DEFAULT_DATASETS_DIR <- "<DATA_ROOT>/realdata/datasets"
RUN_DIR_CANDIDATES <- c("run_CNV", "run_links")
TRUTH_SUBDIRS <- c("wes_cnv", "wgs_cnv")

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

normalize_spot_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("\\.", "-", x)
}

write_standard_csv <- function(mat, out_file) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  write.csv(mat, out_file, quote = FALSE)
}

log_msg <- function(...) {
  cat(sprintf(...), "\n")
}

tool_type_from_name <- function(tool_name) {
  nm <- tolower(tool_name)
  if (grepl("-old$", nm)) return(NA_character_)
  if (nm == "infercnv") return("infercnv")
  if (grepl("^infercnvpy$", nm)) return("infercnvpy")
  if (grepl("^copykat($|_)", nm)) return("copykat")
  if (nm == "scevan") return("scevan")
  if (nm == "starch") return("starch")
  if (nm == "tuscan") return("tuscan")
  if (nm == "conicsmat") return("conicsmat")
  if (nm == "xclone") return("xclone")
  if (nm == "numbat") return("numbat")
  if (nm == "calicost") return("calicost")
  NA_character_
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) return(NULL)
  obj
}

collapse_rows_mean <- function(mat, groups) {
  mat <- as.matrix(mat)
  groups <- as.character(groups)
  sums <- rowsum(mat, group = groups, reorder = FALSE, na.rm = TRUE)
  cnt <- rowsum((!is.na(mat)) * 1, group = groups, reorder = FALSE)
  sums / pmax(cnt, 1)
}

first_matching_name <- function(nms, candidates) {
  nms_key <- normalize_key(nms)
  cand_key <- normalize_key(candidates)
  idx <- match(cand_key, nms_key)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) return(NA_character_)
  nms[idx[1]]
}

copy_number_to_log2 <- function(copy_mat, eps) {
  copy_mat <- as.matrix(copy_mat)
  suppressWarnings(storage.mode(copy_mat) <- "numeric")
  out <- matrix(NA_real_, nrow = nrow(copy_mat), ncol = ncol(copy_mat), dimnames = dimnames(copy_mat))
  ok <- is.finite(copy_mat)
  out[ok] <- log2(pmax(copy_mat[ok], eps) / 2)
  out
}

copy_number_to_state <- function(copy_mat) {
  copy_mat <- as.matrix(copy_mat)
  suppressWarnings(storage.mode(copy_mat) <- "numeric")
  out <- matrix(NA_real_, nrow = nrow(copy_mat), ncol = ncol(copy_mat), dimnames = dimnames(copy_mat))
  ok <- is.finite(copy_mat)
  out[ok] <- ifelse(copy_mat[ok] < 2, -1, ifelse(copy_mat[ok] > 2, 1, 0))
  out
}

state_label_matrix_to_numeric <- function(state_mat) {
  state_chr <- tolower(trimws(as.matrix(state_mat)))
  out <- matrix(0, nrow = nrow(state_chr), ncol = ncol(state_chr), dimnames = dimnames(state_mat))
  out[state_chr %in% c("gain", "amp", "bamp", "hamp", "amplified", "amplification")] <- 1
  out[state_chr %in% c("loss", "del", "bdel", "deletion")] <- -1
  out[state_chr %in% c("neutral", "neu", "loh", "cnloh", "cn_loh")] <- 0

  suppressWarnings(num <- as.numeric(state_chr))
  finite_num <- is.finite(num)
  out[finite_num & num > 0] <- 1
  out[finite_num & num < 0] <- -1
  out[finite_num & num == 0] <- 0
  out
}

normalize_clone_key <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  out <- x_low

  suppressWarnings(num <- as.integer(round(as.numeric(x_low))))
  has_num <- !is.na(num)
  out[has_num] <- paste0("clone", num[has_num])

  prefixed <- grepl("^clone", x_low)
  out[prefixed] <- x_low[prefixed]
  out
}

write_copy_number_outputs <- function(gene_x_spots_copy, tool_dir, eps, prefix) {
  out_log2 <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  out_state <- file.path(tool_dir, "processed_state_spots_x_genes.csv")

  gene_x_spots_copy <- as.matrix(gene_x_spots_copy)
  suppressWarnings(storage.mode(gene_x_spots_copy) <- "numeric")

  if (anyDuplicated(rownames(gene_x_spots_copy)) > 0L) {
    gene_x_spots_copy <- collapse_rows_mean(gene_x_spots_copy, groups = rownames(gene_x_spots_copy))
  }

  final_log2 <- t(copy_number_to_log2(gene_x_spots_copy, eps = eps))
  final_state <- t(copy_number_to_state(gene_x_spots_copy))
  rownames(final_log2) <- normalize_spot_ids(rownames(final_log2))
  rownames(final_state) <- normalize_spot_ids(rownames(final_state))

  write_standard_csv(final_log2, out_log2)
  write_standard_csv(final_state, out_state)
  log_msg("[%s] Saved: %s", prefix, out_log2)
  log_msg("[%s] Saved: %s", prefix, out_state)
  TRUE
}

run_infercnv_pipeline <- function(tool_dir, eps, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_infercnv_expr_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[InferCNV] Skip existing: %s", out_file)
    return(TRUE)
  }

  regular_out <- file.path(paste0(tool_dir, "_regular"), "processed_infercnv_expr_spots_x_genes.csv")
  if (!file.exists(out_file) && file.exists(regular_out)) {
    ok <- file.copy(regular_out, out_file, overwrite = TRUE)
    if (isTRUE(ok)) {
      log_msg("[InferCNV] Copied expr from regular output: %s -> %s", regular_out, out_file)
      return(TRUE)
    }
    log_msg("[InferCNV] Skip: failed to copy regular output %s -> %s", regular_out, out_file)
    return(FALSE)
  }

  f_obj <- file.path(tool_dir, "run.final.infercnv_obj")
  if (!file.exists(f_obj)) {
    # HMM run.final not present -> fall back to the regular run's (identical) expr.data obj.
    regular_obj <- file.path(paste0(tool_dir, "_regular"), "run.final.infercnv_obj")
    if (file.exists(regular_obj)) {
      log_msg("[InferCNV] Using continuous fallback obj: %s", regular_obj)
      f_obj <- regular_obj
    }
  }
  if (!file.exists(f_obj)) {
    log_msg("[InferCNV] Skip: processed log2 not found in %s or _regular, and run.final.infercnv_obj not found in %s or _regular", tool_dir, tool_dir)
    return(FALSE)
  }

  inf_obj <- safe_read_rds(f_obj)
  if (is.null(inf_obj) || inherits(inf_obj, "try-error")) {
    log_msg("[InferCNV] Skip: failed object read in %s", tool_dir)
    return(FALSE)
  }

  pred_mat <- NULL
  if (isS4(inf_obj) && "expr.data" %in% slotNames(inf_obj)) {
    pred_mat <- inf_obj@expr.data
  } else if (is.list(inf_obj) && "expr.data" %in% names(inf_obj)) {
    pred_mat <- inf_obj$expr.data
  }

  if (is.null(pred_mat)) {
    log_msg("[InferCNV] Skip: expr.data not found in %s", tool_dir)
    return(FALSE)
  }

  gene_x_spots <- log2(as.matrix(pred_mat) + eps)
  final_mat <- t(gene_x_spots)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[InferCNV] Saved: %s", out_file)
  TRUE
}

run_infercnv_state_pipeline <- function(tool_dir, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_infercnv_cnv_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[InferCNV-state] Skip existing: %s", out_file)
    return(TRUE)
  }

  candidates <- list.files(tool_dir, pattern = "pred_cnv_genes\\.dat$", full.names = TRUE)
  if (length(candidates) == 0) {
    log_msg("[InferCNV-state] Skip: HMM pred_cnv_genes.dat not found in %s", tool_dir)
    return(FALSE)
  }
  f_state <- sort(candidates)[1]

  dt <- data.table::fread(f_state)
  required_cols <- c("cell_group_name", "gene", "state")
  if (!all(required_cols %in% names(dt))) {
    log_msg("[InferCNV-state] Skip: %s missing %s", f_state, paste(required_cols, collapse = ", "))
    return(FALSE)
  }

  dt <- dt[!is.na(cell_group_name) & !is.na(gene) & !is.na(state)]
  # HMM state -> {-1 loss, 0 neutral, 1 gain}. Mapping depends on HMM_type:
  #   i3: 1=loss, 2=neutral, 3=gain
  #   i6: 1,2=loss, 3=neutral, 4,5,6=gain
  # Detect from inferCNV's filename token (HMMi3 / HMMi6); fall back to max state.
  bn <- basename(f_state)
  if (grepl("HMMi3", bn)) {
    hmm_type <- "i3"
  } else if (grepl("HMMi6", bn)) {
    hmm_type <- "i6"
  } else {
    hmm_type <- if (max(dt$state, na.rm = TRUE) <= 3L) "i3" else "i6"
    log_msg("[InferCNV-state] HMM_type not in filename, inferred %s from max state: %s",
            hmm_type, bn)
  }
  if (hmm_type == "i3") {
    dt[, sval := data.table::fcase(
      state == 1L, -1,
      state == 2L, 0,
      state == 3L, 1,
      default = 0
    )]
  } else {
    dt[, sval := data.table::fcase(
      state %in% c(1L, 2L), -1,
      state == 3L, 0,
      state %in% c(4L, 5L, 6L), 1,
      default = 0
    )]
  }
  log_msg("[InferCNV-state] HMM_type=%s file=%s", hmm_type, bn)
  dt <- dt[, .(sval = sval[1L]), by = .(cell = cell_group_name, gene)]
  mat_w <- data.table::dcast(dt, cell ~ gene, value.var = "sval", fill = 0)
  mat <- as.data.frame(mat_w, check.names = FALSE)
  rownames(mat) <- mat$cell
  mat$cell <- NULL
  mat <- as.matrix(mat)
  rownames(mat) <- normalize_spot_ids(rownames(mat))

  write_standard_csv(mat, out_file)
  log_msg("[InferCNV-state] Saved: %s", out_file)
  TRUE
}

run_infercnvpy_pipeline <- function(tool_dir, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[InferCNVPy] Skip existing: %s", out_file)
    return(TRUE)
  }

  if (!requireNamespace("arrow", quietly = TRUE)) {
    log_msg("[InferCNVPy] Skip: package 'arrow' not available")
    return(FALSE)
  }

  f_parquet <- file.path(tool_dir, "infercnvpy_gene_cnv.parquet")
  if (!file.exists(f_parquet)) {
    log_msg("[InferCNVPy] Skip: infercnvpy_gene_cnv.parquet not found in %s", tool_dir)
    return(FALSE)
  }

  df <- tryCatch(arrow::read_parquet(f_parquet), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) {
    log_msg("[InferCNVPy] Skip: failed read parquet in %s", tool_dir)
    return(FALSE)
  }

  gene_col <- colnames(df)[1]
  mat_df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(mat_df) <- as.character(mat_df[[gene_col]])
  mat_df[[gene_col]] <- NULL
  gene_x_spots <- as.matrix(mat_df)
  suppressWarnings(storage.mode(gene_x_spots) <- "numeric")

  final_mat <- t(gene_x_spots)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[InferCNVPy] Saved: %s", out_file)
  TRUE
}

run_copykat_pipeline <- function(tool_dir, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[CopyKAT] Skip existing: %s", out_file)
    return(TRUE)
  }

  f_raw <- file.path(tool_dir, "_copykat_CNA_raw_results_gene_by_cell.txt")
  if (!file.exists(f_raw)) {
    log_msg("[CopyKAT] Skip: _copykat_CNA_raw_results_gene_by_cell.txt not found in %s", tool_dir)
    return(FALSE)
  }

  dt <- tryCatch(fread(f_raw), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) {
    log_msg("[CopyKAT] Skip: failed read %s", f_raw)
    return(FALSE)
  }

  meta_cols <- c(
    "abspos", "chromosome_name", "start_position", "end_position",
    "ensembl_gene_id", "hgnc_symbol", "band",
    "chrom", "chrompos", "abspos"
  )
  gene_candidates <- c("hgnc_symbol", "gene", "gene_name", "SYMBOL", "symbol", "ensembl_gene_id")
  gene_col <- gene_candidates[gene_candidates %in% names(dt)]
  gene_col <- if (length(gene_col) > 0) gene_col[1] else NA_character_

  if (is.na(gene_col)) {
    log_msg("[CopyKAT] Skip: no gene column found in %s", f_raw)
    return(FALSE)
  }

  genes <- as.character(dt[[gene_col]])
  keep <- !is.na(genes) & nzchar(genes)
  dt <- dt[keep]
  genes <- genes[keep]
  if (length(genes) == 0L) {
    log_msg("[CopyKAT] Skip: no non-empty genes in %s", f_raw)
    return(FALSE)
  }

  cell_cols <- setdiff(names(dt), unique(c(meta_cols, gene_col)))
  if (length(cell_cols) == 0L) {
    log_msg("[CopyKAT] Skip: no cell columns in %s", f_raw)
    return(FALSE)
  }

  gene_x_spots <- as.matrix(dt[, ..cell_cols])
  suppressWarnings(storage.mode(gene_x_spots) <- "numeric")
  rownames(gene_x_spots) <- genes

  if (anyDuplicated(rownames(gene_x_spots)) > 0L) {
    sums <- rowsum(gene_x_spots, group = rownames(gene_x_spots), reorder = FALSE, na.rm = TRUE)
    cnt <- rowsum(!is.na(gene_x_spots), group = rownames(gene_x_spots), reorder = FALSE)
    gene_x_spots <- sums / pmax(cnt, 1)
  }

  final_mat <- t(gene_x_spots)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[CopyKAT] Saved: %s", out_file)
  TRUE
}

run_scevan_pipeline <- function(tool_dir, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[SCEVAN] Skip existing: %s", out_file)
    return(TRUE)
  }

  f_cna <- c(
    file.path(tool_dir, "output", "_CNAmtx.RData"),
    file.path(tool_dir, "_CNAmtx.RData")
  )
  f_cna <- f_cna[file.exists(f_cna)]
  if (length(f_cna) == 0L) {
    log_msg("[SCEVAN] Skip: _CNAmtx.RData not found in %s", tool_dir)
    return(FALSE)
  }

  loaded <- tryCatch(load(f_cna[1]), error = function(e) character(0))
  if (length(loaded) == 0L) {
    log_msg("[SCEVAN] Skip: failed load CNA matrix in %s", tool_dir)
    return(FALSE)
  }
  mat <- get(loaded[1])
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  # Fill gene names from annotation when needed.
  if (is.null(rownames(mat)) || anyNA(rownames(mat)) || any(rownames(mat) == "")) {
    f_annot <- c(
      file.path(tool_dir, "output", "_count_mtx_annot.RData"),
      file.path(tool_dir, "_count_mtx_annot.RData")
    )
    f_annot <- f_annot[file.exists(f_annot)]
    if (length(f_annot) > 0L) {
      loaded_annot <- tryCatch(load(f_annot[1]), error = function(e) character(0))
      if (length(loaded_annot) > 0L) {
        annot <- get(loaded_annot[1])
        genes <- NULL
        if (is.data.frame(annot) && "gene_name" %in% colnames(annot)) genes <- annot$gene_name
        if (is.list(annot) && "gene_name" %in% names(annot)) genes <- annot$gene_name
        if (!is.null(genes) && length(genes) == nrow(mat)) {
          rownames(mat) <- as.character(genes)
        }
      }
    }
  }

  final_mat <- t(mat)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[SCEVAN] Saved: %s", out_file)
  TRUE
}

run_conicsmat_pipeline <- function(tool_dir, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[CONICSmat] Skip existing: %s", out_file)
    return(TRUE)
  }

  f_rds <- file.path(tool_dir, "CONICS_track_matrix.rds")
  if (!file.exists(f_rds)) {
    log_msg("[CONICSmat] Skip: CONICS_track_matrix.rds not found in %s", tool_dir)
    return(FALSE)
  }

  pred_mat <- safe_read_rds(f_rds)
  if (is.null(pred_mat)) {
    log_msg("[CONICSmat] Skip: failed read %s", f_rds)
    return(FALSE)
  }
  if (!is.matrix(pred_mat)) pred_mat <- as.matrix(pred_mat)
  suppressWarnings(storage.mode(pred_mat) <- "numeric")

  final_mat <- t(pred_mat)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[CONICSmat] Saved: %s", out_file)
  TRUE
}

coord_key_from_label <- function(x) {
  parts <- strsplit(as.character(x), "x", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2) return(NA_character_)
    a <- suppressWarnings(as.integer(round(as.numeric(p[1]))))
    b <- suppressWarnings(as.integer(round(as.numeric(p[2]))))
    if (is.na(a) || is.na(b)) return(NA_character_)
    paste0(a, "x", b)
  }, character(1))
}

find_clone_idx <- function(label_keys, state_cols) {
  label_keys <- as.character(label_keys)
  state_cols <- as.character(state_cols)
  if (all(label_keys %in% state_cols)) {
    return(match(label_keys, state_cols))
  }

  lab_num <- suppressWarnings(as.integer(label_keys))
  st_num <- suppressWarnings(as.integer(state_cols))
  if (all(!is.na(lab_num)) && all(!is.na(st_num))) {
    if (all(lab_num %in% st_num)) return(match(lab_num, st_num))
    if (all((lab_num - 1L) %in% st_num)) return(match(lab_num - 1L, st_num))
    if (all((lab_num + 1L) %in% st_num)) return(match(lab_num + 1L, st_num))
  }
  NULL
}

run_starch_pipeline <- function(tool_dir, eps, force = FALSE) {
  out_log2 <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  out_state <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  if (file.exists(out_log2) && file.exists(out_state) && !force) {
    log_msg("[STARCH] Skip existing: %s", tool_dir)
    return(TRUE)
  }

  f_states <- file.path(tool_dir, "states_STARCH.csv")
  f_labels <- file.path(tool_dir, "labels_STARCH.csv")
  if (!file.exists(f_states) || !file.exists(f_labels)) {
    log_msg("[STARCH] Skip: states_STARCH.csv or labels_STARCH.csv missing in %s", tool_dir)
    return(FALSE)
  }

  f_pos <- c(
    file.path(tool_dir, "starch_10x_subset", "tissue_positions_list.csv"),
    file.path(tool_dir, "tissue_positions_list.csv")
  )
  f_pos <- f_pos[file.exists(f_pos)]
  if (length(f_pos) == 0L) {
    log_msg("[STARCH] Skip: tissue_positions_list.csv missing in %s", tool_dir)
    return(FALSE)
  }

  states <- tryCatch(read.csv(f_states, row.names = 1, check.names = FALSE), error = function(e) NULL)
  labels <- tryCatch(read.csv(f_labels, row.names = 1, check.names = FALSE), error = function(e) NULL)
  pos <- tryCatch(fread(f_pos[1], header = FALSE), error = function(e) NULL)
  if (is.null(states) || is.null(labels) || is.null(pos) || ncol(labels) < 1L || ncol(pos) < 4L) {
    log_msg("[STARCH] Skip: failed read states/labels/positions in %s", tool_dir)
    return(FALSE)
  }

  state_mat <- as.matrix(states)
  suppressWarnings(storage.mode(state_mat) <- "numeric")
  clone_labels <- as.character(labels[[1]])
  coord_keys <- coord_key_from_label(rownames(labels))

  # tissue_positions_list.csv: barcode,in_tissue,array_row,array_col,...
  barcodes <- as.character(pos[[1]])
  array_row <- suppressWarnings(as.integer(round(as.numeric(pos[[3]]))))
  array_col <- suppressWarnings(as.integer(round(as.numeric(pos[[4]]))))
  pos_keys <- paste0(array_row, "x", array_col)

  spot_ids <- barcodes[match(coord_keys, pos_keys)]
  clone_idx <- find_clone_idx(clone_labels, colnames(state_mat))

  if (is.null(clone_idx)) {
    log_msg("[STARCH] Skip: failed clone-column mapping in %s", tool_dir)
    return(FALSE)
  }

  keep <- !is.na(spot_ids) & !is.na(clone_idx)
  if (!any(keep)) {
    log_msg("[STARCH] Skip: no valid spot/clone mapping in %s", tool_dir)
    return(FALSE)
  }

  spot_ids <- normalize_spot_ids(spot_ids[keep])
  clone_idx <- clone_idx[keep]

  # STARCH states 0=loss, 1=neutral, 2=gain -> total copy number 1/2/3,
  # then standard log2(cn/2) (-1/0/+0.585) and state (-1/0/+1) via shared helper.
  gene_x_spots_copy <- state_mat[, clone_idx, drop = FALSE] + 1
  colnames(gene_x_spots_copy) <- spot_ids

  write_copy_number_outputs(gene_x_spots_copy, tool_dir = tool_dir, eps = eps, prefix = "STARCH")
}

run_tuscan_pipeline <- function(tool_dir, eps, force = FALSE) {
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  if (file.exists(out_file) && !force) {
    log_msg("[TUSCAN] Skip existing: %s", out_file)
    return(TRUE)
  }

  f_obj <- c(file.path(tool_dir, "cnv_obj.rds.gz"), file.path(tool_dir, "cnv_obj.rds"))
  f_obj <- f_obj[file.exists(f_obj)]
  if (length(f_obj) == 0L) {
    log_msg("[TUSCAN] Skip: cnv_obj.rds.gz/cnv_obj.rds missing in %s", tool_dir)
    return(FALSE)
  }

  obj <- tryCatch({
    if (grepl("\\.gz$", f_obj[1])) readRDS(gzfile(f_obj[1])) else readRDS(f_obj[1])
  }, error = function(e) NULL)
  if (is.null(obj) || inherits(obj, "try-error")) {
    log_msg("[TUSCAN] Skip: failed read %s", f_obj[1])
    return(FALSE)
  }

  cnv_mat <- NULL
  if (isS4(obj) && "cnv.data" %in% slotNames(obj)) {
    cnv_mat <- slot(obj, "cnv.data")
  } else if (is.list(obj) && "cnv.data" %in% names(obj)) {
    cnv_mat <- obj$cnv.data
  }
  if (is.null(cnv_mat)) {
    log_msg("[TUSCAN] Skip: cnv.data not found in %s", f_obj[1])
    return(FALSE)
  }

  gene_x_spots <- as.matrix(cnv_mat)
  suppressWarnings(storage.mode(gene_x_spots) <- "numeric")
  gene_x_spots <- log2(gene_x_spots + eps)

  final_mat <- t(gene_x_spots)
  rownames(final_mat) <- normalize_spot_ids(rownames(final_mat))
  write_standard_csv(final_mat, out_file)
  log_msg("[TUSCAN] Saved: %s", out_file)
  TRUE
}

run_xclone_pipeline <- function(tool_dir, eps, force = FALSE) {
  out_log2 <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  out_state <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  if (file.exists(out_log2) && file.exists(out_state) && !force) {
    log_msg("[XClone] Skip existing: %s", tool_dir)
    return(TRUE)
  }

  f_states <- file.path(tool_dir, "out", "CNV_states_matrix_with_coords_labeled.csv")
  if (!file.exists(f_states)) {
    log_msg("[XClone] Skip: CNV_states_matrix_with_coords_labeled.csv not found in %s", tool_dir)
    return(FALSE)
  }

  dt <- tryCatch(fread(f_states), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L || ncol(dt) < 2L) {
    log_msg("[XClone] Skip: failed read %s", f_states)
    return(FALSE)
  }

  gene_col <- first_matching_name(names(dt), c("GeneName", "gene", "gene_name", "SYMBOL", "symbol"))
  if (is.na(gene_col)) gene_col <- names(dt)[1]

  genes <- as.character(dt[[gene_col]])
  keep_gene <- !is.na(genes) & nzchar(genes)
  if (!any(keep_gene)) {
    log_msg("[XClone] Skip: no non-empty genes in %s", f_states)
    return(FALSE)
  }
  dt <- dt[keep_gene]
  genes <- genes[keep_gene]

  meta_drop <- normalize_key(c(gene_col, "chr", "start", "stop", "geneid", "gene_id"))
  spot_cols <- names(dt)[!(normalize_key(names(dt)) %in% meta_drop)]
  if (length(spot_cols) == 0L) {
    log_msg("[XClone] Skip: no spot columns in %s", f_states)
    return(FALSE)
  }

  state_chr <- as.matrix(dt[, ..spot_cols])
  total_copy <- matrix(
    NA_real_,
    nrow = nrow(state_chr),
    ncol = ncol(state_chr),
    dimnames = list(genes, spot_cols)
  )

  state_low <- tolower(trimws(state_chr))
  total_copy[state_low %in% c("loss", "copyloss")] <- 1
  total_copy[state_low %in% c("neutral", "copyneutral", "loh")] <- 2
  total_copy[state_low %in% c("gain", "copygain")] <- 3

  suppressWarnings(state_num <- as.integer(round(as.numeric(state_chr))))
  total_copy[is.na(total_copy) & !is.na(state_num) & state_num == 0L] <- 1
  total_copy[is.na(total_copy) & !is.na(state_num) & state_num %in% c(1L, 2L)] <- 2
  total_copy[is.na(total_copy) & !is.na(state_num) & state_num == 3L] <- 3

  if (all(is.na(total_copy))) {
    log_msg("[XClone] Skip: could not map XClone states to copy numbers in %s", f_states)
    return(FALSE)
  }

  write_copy_number_outputs(total_copy, tool_dir = tool_dir, eps = eps, prefix = "XClone")
}

# Lazy, cached loader for numbat's bundled hg38 gene annotation (gtf_hg38),
# used for segment->gene mapping (mirrors syntheticdata extract_and_normalize.R).
.get_gtf_hg38 <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    if (!requireNamespace("numbat", quietly = TRUE)) return(NULL)
    e <- new.env(parent = emptyenv())
    ok <- tryCatch({
      utils::data("gtf_hg38", package = "numbat", envir = e)
      TRUE
    }, error = function(err) FALSE)
    if (!ok || !exists("gtf_hg38", envir = e, inherits = FALSE)) return(NULL)
    cache <<- get("gtf_hg38", envir = e, inherits = FALSE)
    cache
  }
})

run_numbat_pipeline <- function(tool_dir, force = FALSE) {
  out_expr <- file.path(tool_dir, "processed_numbat_expr_spots_x_genes.csv")
  out_cnv <- file.path(tool_dir, "processed_numbat_cnv_spots_x_genes.csv")

  normalize_chr <- function(x) {
    out <- tolower(trimws(as.character(x)))
    sub("^chr", "", out)
  }

  # numbat_expr: spot x gene expression signal from gexp_roll_wide.tsv.gz.
  write_numbat_expr <- function() {
    if (file.exists(out_expr) && !force) {
      log_msg("[Numbat-expr] Skip existing: %s", out_expr)
      return(invisible(TRUE))
    }
    f_expr <- file.path(tool_dir, "gexp_roll_wide.tsv.gz")
    if (!file.exists(f_expr)) {
      log_msg("[Numbat-expr] Skip: gexp_roll_wide.tsv.gz missing in %s", tool_dir)
      return(invisible(FALSE))
    }

    dt_expr <- tryCatch(fread(f_expr), error = function(e) NULL)
    if (is.null(dt_expr)) {
      dt_expr <- tryCatch(
        as.data.table(read.delim(
          gzfile(f_expr), sep = "\t", header = TRUE, check.names = FALSE,
          stringsAsFactors = FALSE
        )),
        error = function(e) NULL
      )
    }
    if (is.null(dt_expr) || nrow(dt_expr) == 0L || ncol(dt_expr) < 2L) {
      log_msg("[Numbat-expr] Skip: failed read expression in %s", tool_dir)
      return(invisible(FALSE))
    }
    if (!"cell" %in% names(dt_expr)) setnames(dt_expr, 1L, "cell")

    spot_ids <- normalize_spot_ids(dt_expr$cell)
    mat <- as.matrix(dt_expr[, setdiff(names(dt_expr), "cell"), with = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    rownames(mat) <- spot_ids
    if (anyDuplicated(rownames(mat)) > 0L) {
      mat <- collapse_rows_mean(mat, groups = rownames(mat))
    }

    write_standard_csv(mat, out_expr)
    log_msg("[Numbat-expr] Saved: %s", out_expr)
    invisible(TRUE)
  }

  # numbat_cnv: numeric +1/0/-1 state from joint_post_*.tsv, mapped seg->gene via
  # numbat's bundled gtf_hg38 (mirrors syntheticdata extract_and_normalize.R).
  write_numbat_cnv <- function() {
    if (file.exists(out_cnv) && !force) {
      log_msg("[Numbat-cnv] Skip existing: %s", out_cnv)
      return(invisible(TRUE))
    }

    gtf <- .get_gtf_hg38()
    needed_pkgs <- c("GenomicRanges", "IRanges", "S4Vectors")
    if (is.null(gtf) || !all(vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
      log_msg("[Numbat-cnv] Skip: gtf_hg38/GenomicRanges/IRanges/S4Vectors unavailable in %s", tool_dir)
      return(invisible(FALSE))
    }

    f_joint <- list.files(tool_dir, pattern = "^joint_post_[0-9]+\\.tsv$", full.names = TRUE)
    if (length(f_joint) == 0L) {
      log_msg("[Numbat-cnv] Skip: joint_post_*.tsv missing in %s", tool_dir)
      return(invisible(FALSE))
    }
    joint_iter <- suppressWarnings(as.integer(sub("^joint_post_([0-9]+)\\.tsv$", "\\1", basename(f_joint))))
    f_joint <- f_joint[order(joint_iter, decreasing = TRUE, na.last = TRUE)]

    dt_post <- NULL
    used_joint <- NA_character_
    required_post <- c("cell", "seg", "CHROM", "seg_start", "seg_end", "cnv_state_map")
    for (cand in f_joint) {
      tmp <- tryCatch(fread(cand), error = function(e) NULL)
      if (!is.null(tmp) && all(required_post %in% names(tmp))) {
        dt_post <- tmp
        used_joint <- cand
        break
      }
    }
    if (is.null(dt_post)) {
      log_msg("[Numbat-cnv] Skip: no readable joint_post file with required columns in %s", tool_dir)
      return(invisible(FALSE))
    }

    dt_post <- dt_post[!is.na(cell) & !is.na(seg) & !is.na(CHROM) & !is.na(seg_start) & !is.na(seg_end)]
    dt_post[, `:=`(
      cell = normalize_spot_ids(cell),
      seg = as.character(seg),
      chrom = normalize_chr(CHROM),
      start = suppressWarnings(as.integer(round(as.numeric(seg_start)))),
      end = suppressWarnings(as.integer(round(as.numeric(seg_end)))),
      raw_state = tolower(trimws(as.character(cnv_state_map)))
    )]
    dt_post <- dt_post[!is.na(start) & !is.na(end)]
    if (nrow(dt_post) == 0L) {
      log_msg("[Numbat-cnv] Skip: no usable segments in joint_post in %s", tool_dir)
      return(invisible(FALSE))
    }

    dt_gtf <- data.table::as.data.table(gtf)
    required_gtf <- c("CHROM", "gene_start", "gene_end", "gene")
    if (!all(required_gtf %in% names(dt_gtf))) {
      log_msg("[Numbat-cnv] Skip: gtf_hg38 missing required columns")
      return(invisible(FALSE))
    }
    dt_gtf <- dt_gtf[, .(
      chrom = normalize_chr(CHROM),
      start = suppressWarnings(as.integer(round(as.numeric(gene_start)))),
      end = suppressWarnings(as.integer(round(as.numeric(gene_end)))),
      gene = trimws(as.character(gene))
    )]
    dt_gtf <- dt_gtf[nzchar(gene) & !is.na(start) & !is.na(end)]
    if (nrow(dt_gtf) == 0L) {
      log_msg("[Numbat-cnv] Skip: no usable genes in gtf_hg38")
      return(invisible(FALSE))
    }

    unique_segs <- unique(dt_post[, .(seg, chrom, start, end)])
    gr_segs <- GenomicRanges::GRanges(
      seqnames = unique_segs$chrom,
      ranges = IRanges::IRanges(start = unique_segs$start, end = unique_segs$end),
      seg_id = unique_segs$seg
    )
    gr_genes <- GenomicRanges::GRanges(
      seqnames = dt_gtf$chrom,
      ranges = IRanges::IRanges(start = dt_gtf$start, end = dt_gtf$end),
      gene_name = dt_gtf$gene
    )
    hits <- GenomicRanges::findOverlaps(gr_genes, gr_segs)
    if (length(hits) == 0L) {
      log_msg("[Numbat-cnv] Skip: no segment->gene overlaps against gtf_hg38 in %s", tool_dir)
      return(invisible(FALSE))
    }
    gene_map <- data.table(
      gene = S4Vectors::mcols(gr_genes)$gene_name[S4Vectors::queryHits(hits)],
      seg = S4Vectors::mcols(gr_segs)$seg_id[S4Vectors::subjectHits(hits)]
    )
    gene_map <- unique(gene_map, by = "gene")

    dt <- merge(dt_post[, .(cell, seg, raw_state)], gene_map, by = "seg", allow.cartesian = TRUE)
    if (nrow(dt) == 0L) {
      log_msg("[Numbat-cnv] Skip: joint_post segments did not map to genes via gtf_hg38 in %s", tool_dir)
      return(invisible(FALSE))
    }

    dt[, state_label := data.table::fcase(
      raw_state %in% c("del", "bdel", "loss", "deletion"), "loss",
      raw_state %in% c("loh", "cnloh", "cn_loh"), "loh",
      raw_state %in% c("amp", "bamp", "hamp", "gain", "amplified", "amplification"), "gain",
      raw_state %in% c("neu", "neutral"), "neutral",
      default = "neutral"
    )]
    dt <- unique(dt[, .(cell, gene, state_label)])
    dt <- dt[, .(state_label = state_label[1L]), by = .(cell, gene)]
    mat_w <- data.table::dcast(dt, cell ~ gene, value.var = "state_label", fill = "neutral")
    state_df <- as.data.frame(mat_w, check.names = FALSE)
    rownames(state_df) <- state_df$cell
    state_df$cell <- NULL
    state_mat <- as.matrix(state_df)

    write_standard_csv(state_label_matrix_to_numeric(state_mat), out_cnv)
    log_msg("[Numbat-cnv] Saved numeric state from %s (gtf_hg38 mapping): %s", basename(used_joint), out_cnv)
    invisible(TRUE)
  }

  ok_expr <- write_numbat_expr()
  ok_cnv <- write_numbat_cnv()
  isTRUE(ok_expr) || isTRUE(ok_cnv)
}

run_calicost_pipeline <- function(tool_dir, eps, force = FALSE) {
  out_log2 <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  out_state <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  if (file.exists(out_log2) && file.exists(out_state) && !force) {
    log_msg("[CalicoST] Skip existing: %s", tool_dir)
    return(TRUE)
  }

  f_gene <- file.path(tool_dir, "cnv_genelevel.tsv")
  f_labels <- file.path(tool_dir, "clone_labels.tsv")
  if (!file.exists(f_gene) || !file.exists(f_labels)) {
    log_msg("[CalicoST] Skip: cnv_genelevel.tsv or clone_labels.tsv missing in %s", tool_dir)
    return(FALSE)
  }

  gene_dt <- tryCatch(fread(f_gene), error = function(e) NULL)
  label_dt <- tryCatch(fread(f_labels), error = function(e) NULL)
  if (is.null(gene_dt) || nrow(gene_dt) == 0L || is.null(label_dt) || nrow(label_dt) == 0L) {
    log_msg("[CalicoST] Skip: failed read gene-level or clone-label table in %s", tool_dir)
    return(FALSE)
  }

  gene_col <- first_matching_name(names(gene_dt), c("gene", "GeneName", "gene_name", "SYMBOL", "symbol"))
  if (is.na(gene_col)) gene_col <- names(gene_dt)[1]

  genes <- as.character(gene_dt[[gene_col]])
  keep_gene <- !is.na(genes) & nzchar(genes)
  if (!any(keep_gene)) {
    log_msg("[CalicoST] Skip: no non-empty genes in %s", f_gene)
    return(FALSE)
  }
  gene_dt <- gene_dt[keep_gene]
  genes <- genes[keep_gene]

  orig_cols <- names(gene_dt)
  low_cols <- tolower(orig_cols)
  match_parts <- regmatches(low_cols, regexec("^(clone[^[:space:]]+)\\s+([ab])$", low_cols, perl = TRUE))
  clone_id <- vapply(match_parts, function(x) if (length(x) == 3L) x[2] else NA_character_, character(1))
  allele_id <- vapply(match_parts, function(x) if (length(x) == 3L) toupper(x[3]) else NA_character_, character(1))
  valid_clone_ids <- unique(clone_id[!is.na(clone_id)])
  valid_clone_ids <- valid_clone_ids[nzchar(valid_clone_ids)]

  total_copy_list <- list()
  for (clone_nm in valid_clone_ids) {
    idx_a <- which(clone_id == clone_nm & allele_id == "A")
    idx_b <- which(clone_id == clone_nm & allele_id == "B")
    if (length(idx_a) == 0L || length(idx_b) == 0L) next

    vec_a <- suppressWarnings(as.numeric(gene_dt[[orig_cols[idx_a[1]]]]))
    vec_b <- suppressWarnings(as.numeric(gene_dt[[orig_cols[idx_b[1]]]]))
    total_copy_list[[clone_nm]] <- vec_a + vec_b
  }

  if (length(total_copy_list) == 0L) {
    log_msg("[CalicoST] Skip: no clone A/B column pairs found in %s", f_gene)
    return(FALSE)
  }

  total_copy_by_clone <- do.call(cbind, total_copy_list)
  if (!is.matrix(total_copy_by_clone)) {
    total_copy_by_clone <- matrix(total_copy_by_clone, ncol = 1L)
    colnames(total_copy_by_clone) <- names(total_copy_list)[1]
  }
  rownames(total_copy_by_clone) <- genes

  spot_col <- first_matching_name(names(label_dt), c("BARCODES", "barcode", "barcodes", "spot", "spot_id"))
  clone_col <- first_matching_name(names(label_dt), c("clone_label", "clone", "label"))
  if (is.na(spot_col) || is.na(clone_col)) {
    log_msg("[CalicoST] Skip: missing barcode/clone columns in %s", f_labels)
    return(FALSE)
  }

  spot_ids <- normalize_spot_ids(label_dt[[spot_col]])
  clone_keys <- normalize_clone_key(label_dt[[clone_col]])
  col_idx <- match(clone_keys, colnames(total_copy_by_clone))

  keep <- !is.na(spot_ids) & nzchar(spot_ids) & !is.na(col_idx)
  if (!any(keep)) {
    log_msg("[CalicoST] Skip: no barcode-to-clone matches between %s and %s", f_gene, f_labels)
    return(FALSE)
  }

  spot_ids <- spot_ids[keep]
  col_idx <- col_idx[keep]
  gene_x_spots_copy <- total_copy_by_clone[, col_idx, drop = FALSE]
  colnames(gene_x_spots_copy) <- spot_ids

  if (anyDuplicated(colnames(gene_x_spots_copy)) > 0L) {
    spots_x_genes_copy <- t(gene_x_spots_copy)
    spots_x_genes_copy <- collapse_rows_mean(spots_x_genes_copy, groups = rownames(spots_x_genes_copy))
    gene_x_spots_copy <- t(spots_x_genes_copy)
  }

  write_copy_number_outputs(gene_x_spots_copy, tool_dir = tool_dir, eps = eps, prefix = "CalicoST")
}

process_one_tool <- function(tool_dir, eps, force = FALSE) {
  tool_name <- basename(tool_dir)
  tool_type <- tool_type_from_name(tool_name)
  if (is.na(tool_type)) return(invisible(NULL))

  if (tool_type == "infercnv") {
    run_infercnv_pipeline(tool_dir, eps = eps, force = force)
    run_infercnv_state_pipeline(tool_dir, force = force)
  } else if (tool_type == "infercnvpy") {
    run_infercnvpy_pipeline(tool_dir, force = force)
  } else if (tool_type == "copykat") {
    run_copykat_pipeline(tool_dir, force = force)
  } else if (tool_type == "scevan") {
    run_scevan_pipeline(tool_dir, force = force)
  } else if (tool_type == "starch") {
    run_starch_pipeline(tool_dir, eps = eps, force = force)
  } else if (tool_type == "tuscan") {
    run_tuscan_pipeline(tool_dir, eps = eps, force = force)
  } else if (tool_type == "conicsmat") {
    run_conicsmat_pipeline(tool_dir, force = force)
  } else if (tool_type == "xclone") {
    run_xclone_pipeline(tool_dir, eps = eps, force = force)
  } else if (tool_type == "numbat") {
    run_numbat_pipeline(tool_dir, force = force)
  } else if (tool_type == "calicost") {
    run_calicost_pipeline(tool_dir, eps = eps, force = force)
  }
}

find_run_dir <- function(dataset_dir) {
  for (cand in RUN_DIR_CANDIDATES) {
    p <- file.path(dataset_dir, cand)
    if (dir.exists(p)) return(p)
  }
  NA_character_
}

discover_dataset_roots <- function(datasets_dir, datasets_filter = NULL) {
  ds_dirs <- list.dirs(datasets_dir, recursive = FALSE, full.names = TRUE)
  ds_dirs <- ds_dirs[file.info(ds_dirs)$isdir]
  if (!is.null(datasets_filter)) {
    ds_dirs <- ds_dirs[normalize_key(basename(ds_dirs)) %in% normalize_key(datasets_filter)]
  }
  out <- list()
  for (d in ds_dirs) {
    truth_dir <- file.path(d, "truth_links")
    has_truth <- any(vapply(TRUTH_SUBDIRS, function(s) dir.exists(file.path(truth_dir, s)), logical(1)))
    if (!has_truth) next
    run_dir <- find_run_dir(d)
    if (is.na(run_dir)) {
      log_msg("[Init] Skip dataset %s: no %s present", basename(d), paste(RUN_DIR_CANDIDATES, collapse = "/"))
      next
    }
    out[[basename(d)]] <- run_dir
  }
  out
}

process_root <- function(root_dir, dataset_name, samples, strategies, tools, eps, force) {
  sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[file.info(sample_dirs)$isdir]
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[normalize_key(basename(sample_dirs)) %in% normalize_key(samples)]
  }

  log_msg("\n[Dataset] %s | root_dir = %s | samples = %d",
          dataset_name, root_dir, length(sample_dirs))

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    run_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[file.info(run_dirs)$isdir]
    if (!is.null(strategies)) {
      run_dirs <- run_dirs[normalize_key(basename(run_dirs)) %in% normalize_key(strategies)]
    }
    log_msg("\n[Sample] %s | runs = %d", sample_name, length(run_dirs))

    for (run_dir in run_dirs) {
      run_name <- basename(run_dir)
      tool_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
      tool_dirs <- tool_dirs[file.info(tool_dirs)$isdir]
      tool_dirs <- tool_dirs[!tolower(basename(tool_dirs)) %in% c("figures", "tables")]
      log_msg("  [Run] %s | tools(found dirs) = %d", run_name, length(tool_dirs))

      for (tool_dir in tool_dirs) {
        tool_name <- basename(tool_dir)
        tool_type <- tool_type_from_name(tool_name)
        if (is.na(tool_type)) next
        if (!is.null(tools)) {
          tool_keys <- normalize_key(tools)
          if (!(normalize_key(tool_name) %in% tool_keys || normalize_key(tool_type) %in% tool_keys)) {
            next
          }
        }
        log_msg("    [Tool] %s (%s)", tool_name, tool_type)

        tryCatch(
          process_one_tool(tool_dir, eps = eps, force = force),
          error = function(e) {
            log_msg("    [Tool] ERROR in %s: %s", tool_name, e$message)
          }
        )
      }
    }
  }
}

main <- function() {
  datasets_dir <- DEFAULT_DATASETS_DIR
  root_dir <- NULL
  datasets <- NULL
  samples <- NULL
  strategies <- NULL
  tools <- NULL
  force <- FALSE
  eps <- EPS

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L) {
    i <- 1L
    while (i <= length(args)) {
      key <- args[i]
      if (key == "--datasets_dir" && i + 1L <= length(args)) {
        datasets_dir <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--root_dir" && i + 1L <= length(args)) {
        root_dir <- args[i + 1L]; i <- i + 2L; next
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
      if (key == "--eps" && i + 1L <= length(args)) {
        eps <- as.numeric(args[i + 1L]); i <- i + 2L; next
      }
      if (key == "--force") {
        force <- TRUE; i <- i + 1L; next
      }
      i <- i + 1L
    }
  }

  if (!is.finite(eps) || eps <= 0) stop("--eps must be > 0")

  if (!is.null(root_dir)) {
    if (!dir.exists(root_dir)) stop("root_dir does not exist: ", root_dir)
    roots <- setNames(list(root_dir), basename(dirname(normalizePath(root_dir))))
  } else {
    if (!dir.exists(datasets_dir)) stop("datasets_dir does not exist: ", datasets_dir)
    roots <- discover_dataset_roots(datasets_dir, datasets_filter = datasets)
    if (length(roots) == 0L) stop("No datasets with WES/WGS truth found under: ", datasets_dir)
  }

  log_msg("[Init] datasets = %d (%s)", length(roots), paste(names(roots), collapse = ", "))
  if (!is.null(strategies)) log_msg("[Init] strategies filter = %s", paste(strategies, collapse = ", "))
  if (!is.null(tools)) log_msg("[Init] tools filter = %s", paste(tools, collapse = ", "))
  log_msg("[Init] force = %s", as.character(force))

  for (ds_name in names(roots)) {
    process_root(
      root_dir = roots[[ds_name]], dataset_name = ds_name,
      samples = samples, strategies = strategies, tools = tools,
      eps = eps, force = force
    )
  }

  log_msg("\n[Done] Standardization finished.")
}

main()
