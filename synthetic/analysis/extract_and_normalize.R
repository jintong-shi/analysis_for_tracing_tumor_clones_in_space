#!/usr/bin/env Rscript
# Extract each tool's raw output from Tool_Outputs and normalize it to a common
# Spots x Genes log2 (or +1/0/-1 state) CSV, one per tool per Rep_* batch.
library(Matrix)
library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)
library(dplyr)
library(arrow)
library(tibble)
EPS=1e-12
library(EnsDb.Hsapiens.v86)
library(numbat)
data("gtf_hg38", package = "numbat")
HG38_CHROMINFO_FILE <- "<REF_ROOT>/hg38_chromInfo.txt.gz"

# Which output root to process under each batch, e.g. "Tool_Outputs" or "Tool_Output_TME"
TOOL_OUT_NAME <- "Tool_Outputs"
# Suffix appended to tool folders, e.g. "_TME" for CopyKAT_TME / Numbat_TME.
TOOL_OUT_SUFFIX <- ""
# InferCNV folder name under TOOL_OUT_NAME, e.g. "InferCNV" or "InferCNV_noHMM"
INFERCNV_DIR_NAME <- "InferCNV"
get_tool_dir <- function(batch_path, tool_name) {
  file.path(batch_path, TOOL_OUT_NAME, paste0(tool_name, TOOL_OUT_SUFFIX))
}

# ==============================================================================
# Helper: state label -> numeric matrix
# gains -> +1 ; neu/loh -> 0 ; loss -> -1
# ==============================================================================
state_matrix_to_numeric <- function(state_mat) {
  x_df <- as.data.frame(state_mat, stringsAsFactors = FALSE, check.names = FALSE)
  row_ids <- rownames(x_df)
  
  to_num <- function(x) {
    s <- tolower(trimws(as.character(x)))
    out <- rep(0, length(s))
    
    # Explicit label mapping
    out[s %in% c("gain", "amp", "bamp", "hamp", "amplified", "amplification")] <- 1
    out[s %in% c("loss", "del", "bdel", "deletion")] <- -1
    out[s %in% c("neutral", "neu", "loh", "cnloh", "cn_loh")] <- 0
    
    # Numeric fallback (if any tool writes numeric-like states)
    v <- suppressWarnings(as.numeric(s))
    is_num <- !is.na(v)
    out[is_num & v > 0] <- 1
    out[is_num & v < 0] <- -1
    out[is_num & v == 0] <- 0
    out
  }
  
  num_df <- as.data.frame(lapply(x_df, to_num), check.names = FALSE)
  rownames(num_df) <- row_ids
  as.matrix(num_df)
}

# ==============================================================================
# InferCNV
# ==============================================================================
run_infercnv_pipeline <- function(batch_path,EPS) {
  batch_name <- basename(batch_path)
  message(sprintf("[InferCNV] Processing %s ...", batch_name))
  
  # Paths
  tool_dir <- get_tool_dir(batch_path, INFERCNV_DIR_NAME)
  
  if (!dir.exists(tool_dir)) return(NULL)
  
  f_obj <- file.path(tool_dir, 'run.final.infercnv_obj')
  # HMM run.final not ready yet -> fall back to the regular run's (identical) expr.data
  if (!file.exists(f_obj)) {
    fallback_tool_names <- paste0(INFERCNV_DIR_NAME, "_regular")
    if (grepl("HMM", INFERCNV_DIR_NAME)) {
      fallback_tool_names <- c(
        fallback_tool_names,
        sub("_HMM[0-9]+", "_noHMM", INFERCNV_DIR_NAME)
      )
    }
    for (fallback_tool in unique(fallback_tool_names)) {
      fallback_obj <- file.path(get_tool_dir(batch_path, fallback_tool), 'run.final.infercnv_obj')
      if (file.exists(fallback_obj)) {
        message("  [InferCNV] Using continuous fallback: ", fallback_obj)
        f_obj <- fallback_obj
        break
      }
    }
  }
  f_txt <- file.path(tool_dir, 'infercnv.observations.txt')
  
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  message(sprintf("  Processing InferCNV in: %s", basename(batch_path)))
  
  # Read the object if present, otherwise the txt
  if (!file.exists(f_obj)) stop('inferCNV prediction file not found: ', f_obj)
  inf_obj <- readRDS(f_obj)
  pred_mat <- inf_obj@expr.data
  
  # Log2
  gene_x_cells <- log2(pred_mat + EPS)
  
  # Transpose to Spots x Genes and save
  final_mat <- t(gene_x_cells)
  rownames(final_mat) <- sub("^Rep_[0-9]+__", "", rownames(final_mat))
  write.csv(final_mat, out_file, quote = FALSE)
  
  message("[InferCNV] Done. Saved to: ", out_file)
}

run_infercnv_state_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[InferCNV-state] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, INFERCNV_DIR_NAME)
  if (!dir.exists(tool_dir)) return(NULL)
  
  candidates <- list.files(
    tool_dir,
    pattern = "pred_cnv_genes\\.dat$",
    full.names = TRUE
  )
  if (length(candidates) == 0) {
    message("  Skip: InferCNV HMM state file not found.")
    return(NULL)
  }
  f_state <- sort(candidates)[1]
  
  out_file <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  
  dt <- data.table::fread(f_state)
  required_cols <- c("cell_group_name", "gene", "state")
  if (!all(required_cols %in% names(dt))) {
    stop("InferCNV HMM state file must contain: ", paste(required_cols, collapse = ", "))
  }
  
  dt <- dt[!is.na(cell_group_name) & !is.na(gene) & !is.na(state)]
  # HMM state -> {loss, neutral, gain}. Mapping depends on HMM_type:
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
    message(sprintf("  [InferCNV-state] HMM_type not in filename, inferred %s from max state: %s",
                    hmm_type, bn))
  }
  if (hmm_type == "i3") {
    dt[, label := data.table::fcase(
      state == 1L, "loss",
      state == 2L, "neutral",
      state == 3L, "gain",
      default = "neutral"
    )]
  } else {
    dt[, label := data.table::fcase(
      state %in% c(1L, 2L), "loss",
      state == 3L, "neutral",
      state %in% c(4L, 5L, 6L), "gain",
      default = "neutral"
    )]
  }
  message(sprintf("  [InferCNV-state] HMM_type=%s file=%s", hmm_type, bn))
  
  dt <- dt[, .(label = label[1L]), by = .(cell = cell_group_name, gene)]
  mat_labels <- data.table::dcast(dt, cell ~ gene, value.var = "label", fill = "neutral")
  mat <- as.data.frame(mat_labels, check.names = FALSE)
  rownames(mat) <- mat$cell
  mat$cell <- NULL
  
  write.csv(mat, out_file, quote = FALSE)
  message("[InferCNV-state] Done. Saved to: ", out_file)
}

run_copykat_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[CopyKAT] Processing %s ...", batch_name))

  # A. Define Paths
  tool_dir <- get_tool_dir(batch_path, "CopyKAT")
  
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  f_rds    <- file.path(tool_dir, 'copykat_object.rds')
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  # B. Load Data
  if (!file.exists(f_rds)) {
    # If the RDS is missing, we can't proceed (unlike InferCNV where we had txt fallback)
    message(sprintf("  Skip: CopyKAT result object not found in %s", tool_dir))
    return(NULL)
  }
  
  res <- readRDS(f_rds)
  
  # Check if CNAmat exists (it might be missing if CopyKAT failed to find aneuploidy)
  if (is.null(res$CNAmat)) {
    message("  Skip: 'CNAmat' missing in CopyKAT results (possibly no aneuploidy detected).")
    return(NULL)
  }
  
  # C. Extract Matrix
  # CopyKAT CNAmat structure: Rows=Genes, Cols=Info(3cols) + Cells
  # Columns 1-3 are usually: chrom, chrompos, abspos. We remove them.
  # The values are already Log2 Ratios (neutral = 0).
  raw_mat <- res$CNAmat
  
  # Position -> GRanges
  pos <- raw_mat[, c('chrom','chrompos')]
  pos$chrom <- paste0('chr', pos$chrom)
  pos$chrom <- gsub('chr23','chrX', pos$chrom)

  if (!file.exists(HG38_CHROMINFO_FILE)) {
    stop("Local hg38 chromInfo file not found: ", HG38_CHROMINFO_FILE)
  }
  chrom_info <- data.table::fread(HG38_CHROMINFO_FILE, header = FALSE)
  standard_chr <- intersect(
    chrom_info[[1]],
    c(paste0("chr", 1:22), "chrX", "chrY")
  )
  
  edb <- EnsDb.Hsapiens.v86
  genes_gr_global <- genes(edb)
  # Avoid GenomeInfoDb's online UCSC lookup; EnsDb seqlevels are already hg38.
  renamed_seqlevels <- seqlevels(genes_gr_global)
  renamed_seqlevels <- ifelse(grepl("^chr", renamed_seqlevels),
                              renamed_seqlevels,
                              paste0("chr", renamed_seqlevels))
  renamed_seqlevels[renamed_seqlevels == "chrMT"] <- "chrM"
  seqlevels(genes_gr_global) <- renamed_seqlevels
  
  genes_gr_global <- keepSeqlevels(genes_gr_global, standard_chr, pruning.mode = "coarse")
  mcols(genes_gr_global)$gene <- mcols(genes_gr_global)$gene_name
  
  gr_cnv <- GRanges(seqnames=pos$chrom,
                    ranges=IRanges(start=pos$chrompos, end=pos$chrompos))
  # Map to the global gene annotation
  gr_genes <- keepSeqlevels(genes_gr_global, seqlevels(gr_cnv), pruning.mode='coarse')
  common   <- intersect(seqlevels(gr_cnv), seqlevels(gr_genes))
  gr_cnv   <- keepSeqlevels(gr_cnv, common, pruning.mode='coarse')
  gr_genes <- keepSeqlevels(gr_genes, common, pruning.mode='coarse')
  hits <- findOverlaps(gr_cnv, gr_genes)
  pos$gene <- NA
  pos$gene[queryHits(hits)] <- mcols(gr_genes)$gene[subjectHits(hits)]
  mat_full <- cbind(pos, raw_mat[ , -(1:3)])
  mat_filt <- mat_full %>% dplyr::filter(!is.na(gene)) %>% dplyr::distinct(gene, .keep_all=TRUE)
  gene_x_cells <- mat_filt %>% tibble::column_to_rownames('gene') %>% dplyr::select(starts_with('Spot_')) %>% as.matrix()
  
  # D. Transpose to Spots x Genes and Save
  final_mat <- t(as.matrix(gene_x_cells))
  
  write.csv(final_mat, out_file, quote = FALSE)
  
  message("[CopyKAT] Done. Saved to: ", out_file)
}

run_infercnvpy_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[InferCNVPy] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "InferCNVPy")
  if (!dir.exists(tool_dir)) return(NULL)
  
  # InferCNVPy output is usually named "infercnvpy_gene_cnv.parquet"
  f_parquet <- file.path(tool_dir, "infercnvpy_gene_cnv.parquet")
  out_file  <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  if (!file.exists(f_parquet)) {
    message("  Skip: Result parquet not found.")
    return(NULL)
  }
  
  # Read Parquet
  # The Python script outputs: index=gene, cols=cells.
  df <- arrow::read_parquet(f_parquet)
  
  # The first column is usually the gene name (because index was reset in Python)
  gene_col <- colnames(df)[1] # likely "gene"
  
  # Convert to matrix
  mat <- as.data.frame(df)
  rownames(mat) <- mat[[gene_col]]
  mat[[gene_col]] <- NULL
  final_mat <- as.matrix(mat)
  
  # Check format: We need Rows=Spots, Cols=Genes
  # Currently final_mat is Rows=Genes, Cols=Spots
  final_mat_t <- t(final_mat)
  
  # Write CSV
  write.csv(final_mat_t, out_file, quote = FALSE)
  message("[InferCNVPy] Done. Saved to: ", out_file)
}

# ==============================================================================
# SCEVAN
# ==============================================================================
run_scevan_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[SCEVAN] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "SCEVAN")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  f_cna   <- file.path(tool_dir,"output", "_CNAmtx.RData")
  f_annot <- file.path(tool_dir,"output", "_count_mtx_annot.RData")
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  if (!file.exists(f_cna) || !file.exists(f_annot)) {
    message("  Skip: SCEVAN RData files not found.")
    return(NULL)
  }
  
  cnv   <- get(load(f_cna))           # CNA_mtx_relat
  annot <- get(load(f_annot))         # count_mtx_annot
  mat <- cnv
  if (!is.null(annot$gene_name) && length(annot$gene_name) == nrow(mat)) {
    rownames(mat) <- annot$gene_name
  } else {
    warning("SCEVAN: gene_name not found or length mismatch; rownames not set.")
  }
  
  final_mat <- t(mat)  # Spots x Genes
  write.csv(final_mat, out_file, quote = FALSE)
  message("[SCEVAN] Done. Saved to: ", out_file)
}

# ==============================================================================
# CONICSmat
# ==============================================================================
run_conicsmat_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[CONICSmat] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "CONICSmat")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  f_rds <- file.path(tool_dir, "CONICS_track_matrix.rds")
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  if (!file.exists(f_rds)) {
    message("  Skip: CONICSmat result object not found.")
    return(NULL)
  }
  
  pred_mat <- readRDS(f_rds)
  if (!is.matrix(pred_mat)) pred_mat <- as.matrix(pred_mat)
  
  # CONICSmat output: Genes x Cells -> transpose to Spots x Genes
  final_mat <- t(pred_mat)
  write.csv(final_mat, out_file, quote = FALSE)
  message("[CONICSmat] Done. Saved to: ", out_file)
}

# ==============================================================================
# STARCH
# ==============================================================================
run_starch_pipeline <- function(batch_path, EPS=1e-12) {
  batch_name <- basename(batch_path)
  message(sprintf("[STARCH] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "STARCH")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  f_states <- file.path(tool_dir, "states_STITCH_output.csv")
  f_labels <- file.path(tool_dir, "labels_STITCH_output.csv")
  out_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  if (!file.exists(f_states) || !file.exists(f_labels)) {
    message("  Skip: STARCH output files not found.")
    return(NULL)
  }
  
  states <- read.csv(f_states, row.names = 1, check.names = FALSE)
  labels_df <- read.csv(f_labels, row.names = 1, check.names = FALSE)
  if (ncol(labels_df) < 1) {
    message("  Skip: labels file has no columns.")
    return(NULL)
  }
  
  labels_vec <- labels_df[[1]]
  label_coords <- rownames(labels_df)
  
  # Load spatial metadata to map "x x y" back to Spot_* IDs
  meta_path <- file.path(batch_path, "synthetic", "spatial_data", "st-meta.parquet")
  if (!file.exists(meta_path)) {
    message("  Skip: st-meta.parquet not found (cannot map spot labels).")
    return(NULL)
  }
  meta <- as.data.frame(arrow::read_parquet(meta_path))
  if (!all(c("index","x","y") %in% names(meta))) {
    stop("st-meta.parquet must contain columns: index, x, y")
  }
  
  # Replicate prep_starch.py coordinate naming: int(round(x))xint(round(y))
  coord <- paste0(as.integer(round(meta$x)), "x", as.integer(round(meta$y)))
  
  # Handle duplicated coordinate names the same way as in prep_starch.py
  if (any(duplicated(coord))) {
    seen <- new.env(parent = emptyenv())
    new_names <- character(length(coord))
    for (i in seq_along(coord)) {
      name <- coord[i]
      if (!exists(name, envir = seen, inherits = FALSE)) {
        assign(name, 0, envir = seen)
        new_names[i] <- name
      } else {
        n <- get(name, envir = seen, inherits = FALSE) + 1
        assign(name, n, envir = seen)
        new_names[i] <- paste0(name, "_", n)
      }
    }
    coord <- new_names
  }
  
  # Normalize possible "2.0x2.0" formatting from STARCH outputs
  normalize_coord <- function(x) gsub("\\.0", "", x)
  coord_norm <- normalize_coord(coord)
  label_coords_norm <- normalize_coord(label_coords)
  
  meta_map <- data.frame(coord = coord_norm, spot_id = meta$index, stringsAsFactors = FALSE)
  spot_ids <- meta_map$spot_id[match(label_coords_norm, meta_map$coord)]
  
  if (any(is.na(spot_ids))) {
    warning(sum(is.na(spot_ids)), " STARCH labels could not be mapped to Spot_* IDs; dropping.")
  }
  
  states_mat <- as.matrix(states)
  col_names <- colnames(states_mat)
  label_keys <- as.character(labels_vec)
  
  # Determine clone-column indices
  clone_idx <- NULL
  if (all(label_keys %in% col_names)) {
    clone_idx <- match(label_keys, col_names)
  } else {
    col_nums <- suppressWarnings(as.integer(col_names))
    lab_nums <- suppressWarnings(as.integer(label_keys))
    if (!all(is.na(col_nums)) && !all(is.na(lab_nums))) {
      if (all(lab_nums %in% col_nums)) {
        clone_idx <- match(lab_nums, col_nums)
      } else if (all((lab_nums - 1) %in% col_nums)) {
        clone_idx <- match(lab_nums - 1, col_nums)
      }
    }
  }
  
  if (is.null(clone_idx) || any(is.na(clone_idx))) {
    stop("Cannot map STARCH labels to clone columns in states file.")
  }
  
  valid <- !is.na(spot_ids) & !is.na(clone_idx)
  if (!all(valid)) {
    spot_ids <- spot_ids[valid]
    clone_idx <- clone_idx[valid]
  }
  
  # Broadcast clone-level states to spot-level
  gene_x_spots <- states_mat[, clone_idx, drop = FALSE]
  colnames(gene_x_spots) <- spot_ids
  
  # STARCH states 0=loss, 1=neutral, 2=gain -> total copy number 1/2/3.
  total_copy <- gene_x_spots + 1

  # Standard log2(cn/2): loss -1 / neutral 0 / gain +0.585
  log2_mat <- log2(total_copy / 2)
  final_mat <- t(log2_mat)  # Spots x Genes
  write.csv(final_mat, out_file, quote = FALSE)

  # Standard discrete state -1/0/+1 (same convention as other tools)
  state_std <- ifelse(total_copy < 2, -1, ifelse(total_copy > 2, 1, 0))
  out_state_file <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  write.csv(t(state_std), out_state_file, quote = FALSE)

  message("[STARCH] Done. log2 -> ", out_file)
  message("[STARCH] Done. state (-1/0/+1) -> ", out_state_file)
}

# ==============================================================================
# XClone (state matrix)
# ==============================================================================
run_xclone_state_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[XClone] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "XClone")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  # Prefer prob1_merge if present, otherwise any *_cnv_state_spot_by_gene.tsv
  f_preferred <- file.path(tool_dir, "prob1_merge_cnv_state_spot_by_gene.tsv")
  if (file.exists(f_preferred)) {
    f_state <- f_preferred
  } else {
    candidates <- list.files(tool_dir, pattern = "_cnv_state_spot_by_gene\\.tsv$", full.names = TRUE)
    if (length(candidates) == 0) {
      message("  Skip: CNV state TSV not found.")
      return(NULL)
    }
    f_state <- candidates[1]
  }
  
  out_state_file <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  out_num_file   <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  
  # TSV is already Spots x Genes (rows=spots, cols=genes)
  mat <- read.delim(f_state, row.names = 1, check.names = FALSE)
  write.csv(mat, out_state_file, quote = FALSE)
  
  # Also save numeric version for metric compatibility
  num_mat <- state_matrix_to_numeric(mat)
  write.csv(num_mat, out_num_file, quote = FALSE)
  
  message("[XClone] Done. Saved to: ", out_state_file)
  message("[XClone] Numeric state (+1/0/-1) saved to: ", out_num_file)
}

# ==============================================================================
# Numbat (state matrix, Spots x Genes)
# ==============================================================================
run_numbat_state_pipeline <- function(batch_path) {
  batch_name <- basename(batch_path)
  message(sprintf("[Numbat] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "Numbat")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  out_state_file <- file.path(tool_dir, "processed_state_spots_x_genes.csv")
  out_num_file   <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  f_expr <- file.path(tool_dir, "gexp_roll_wide.tsv.gz")
  
  write_numbat_log2_fallback <- function() {
    if (!file.exists(f_expr)) {
      stop("Numbat fallback failed: gexp_roll_wide.tsv.gz not found.")
    }
    dt_expr <- data.table::fread(f_expr)
    if (ncol(dt_expr) < 2) {
      stop("Numbat fallback failed: gexp_roll_wide.tsv.gz must contain cell column plus gene columns.")
    }
    if (!"cell" %in% names(dt_expr)) {
      data.table::setnames(dt_expr, 1L, "cell")
    }
    mat_expr <- as.data.frame(dt_expr, check.names = FALSE)
    rownames(mat_expr) <- mat_expr$cell
    mat_expr$cell <- NULL
    write.csv(mat_expr, out_num_file, quote = FALSE)
    message("[Numbat] Fallback log2 matrix saved from gexp_roll_wide.tsv.gz to: ", out_num_file)
    invisible(mat_expr)
  }
  
  tryCatch({
    if (!requireNamespace("numbat", quietly = TRUE)) {
      stop("Numbat package not available. Please install/load 'numbat'.")
    }
    
    joint_files <- list.files(
      tool_dir,
      pattern = "^joint_post_[0-9]+\\.tsv$",
      full.names = TRUE
    )
    if (length(joint_files) == 0) {
      stop("No Numbat joint_post_*.tsv files found.")
    }
    joint_ids <- as.integer(sub("^joint_post_([0-9]+)\\.tsv$", "\\1", basename(joint_files)))
    joint_files <- joint_files[order(joint_ids, decreasing = TRUE)]
    required_post <- c("cell", "seg", "CHROM", "seg_start", "seg_end", "cnv_state_map")
    read_errors <- character()
    dt_post <- NULL
    chosen_joint <- NULL
    for (joint_file in joint_files) {
      dt_try <- tryCatch(
        data.table::fread(joint_file),
        error = function(e) {
          read_errors <<- c(read_errors, sprintf("%s: %s", basename(joint_file), e$message))
          NULL
        }
      )
      if (is.null(dt_try)) next
      missing_cols <- setdiff(required_post, names(dt_try))
      if (length(missing_cols) > 0) {
        read_errors <- c(
          read_errors,
          sprintf("%s missing columns: %s", basename(joint_file), paste(missing_cols, collapse = ", "))
        )
        next
      }
      dt_post <- dt_try
      chosen_joint <- basename(joint_file)
      break
    }
    if (is.null(dt_post)) {
      stop("No readable Numbat joint_post_*.tsv with required columns. ",
           paste(read_errors, collapse = "; "))
    }
    message("[Numbat] Using ", chosen_joint, " for state extraction.")
    dt_post <- dt_post[, .(cell, seg, chrom = CHROM, start = seg_start, end = seg_end, raw_state = cnv_state_map)]
    
    dt_gtf <- data.table::as.data.table(gtf_hg38)
    required_gtf <- c("CHROM", "gene_start", "gene_end", "gene")
    missing_gtf <- setdiff(required_gtf, names(dt_gtf))
    if (length(missing_gtf) > 0) {
      stop("gtf_hg38 missing columns: ", paste(missing_gtf, collapse = ", "))
    }
    dt_gtf <- dt_gtf[, .(chrom = CHROM, start = gene_start, end = gene_end, gene)]
    
    # Normalize chromosome naming
    norm_chr <- function(x) gsub("^chr", "", as.character(x))
    dt_post[, chrom := norm_chr(chrom)]
    dt_gtf[, chrom := norm_chr(chrom)]
    
    # Unique segments
    unique_segs <- unique(dt_post[, .(seg, chrom, start, end)])
    
    gr_segs <- GRanges(
      seqnames = unique_segs$chrom,
      ranges = IRanges(start = unique_segs$start, end = unique_segs$end),
      seg_id = unique_segs$seg
    )
    
    gr_genes <- GRanges(
      seqnames = dt_gtf$chrom,
      ranges = IRanges(start = dt_gtf$start, end = dt_gtf$end),
      gene_name = dt_gtf$gene
    )
    hits <- findOverlaps(gr_genes, gr_segs)
    gene_map <- data.table(
      gene = mcols(gr_genes)$gene_name[queryHits(hits)],
      seg = mcols(gr_segs)$seg_id[subjectHits(hits)]
    )
    gene_map <- unique(gene_map, by = "gene")
    
    final_dt <- merge(dt_post[, .(cell, seg, raw_state)], gene_map, by = "seg", allow.cartesian = TRUE)
    if (nrow(final_dt) == 0) {
      stop("No segment->gene mapping after merge.")
    }
    
    # Map raw state to normalized labels
    final_dt[, label := data.table::fcase(
      raw_state %in% c("del", "bdel"), "loss",
      raw_state == "loh", "loh",
      raw_state == "neu", "neutral",
      raw_state %in% c("amp", "bamp", "hamp"), "gain",
      default = "neutral"
    )]
    
    mat_labels <- data.table::dcast(final_dt, cell ~ gene, value.var = "label", fill = "neutral")
    mat <- as.data.frame(mat_labels)
    rownames(mat) <- mat$cell
    mat$cell <- NULL
    
    write.csv(mat, out_state_file, quote = FALSE)
    
    # Prefer state-derived numeric labels when available.
    num_mat <- state_matrix_to_numeric(mat)
    write.csv(num_mat, out_num_file, quote = FALSE)
    
    message("[Numbat] Done. Saved to: ", out_state_file)
    message("[Numbat] Numeric state (+1/0/-1) saved to: ", out_num_file)
  }, error = function(e) {
    message("[Numbat] State extraction failed: ", e$message)
    write_numbat_log2_fallback()
  })
}

# ==============================================================================
# CalicoST
# - Sum clone A/B to total copy number
# - Broadcast clone-level CN to spots using clone_labels.tsv
# - Save requested gene x spot total-CN matrix and clone labels
# - Keep legacy processed_log2_spots_x_genes.csv for downstream evaluation
# ==============================================================================
run_calicost_pipeline <- function(batch_path, EPS=1e-12) {
  batch_name <- basename(batch_path)
  message(sprintf("[CalicoST] Processing %s ...", batch_name))
  
  tool_dir <- get_tool_dir(batch_path, "CalicoST")
  if (!dir.exists(tool_dir)) {
    message("  Skip: Output directory not found.")
    return(NULL)
  }
  
  # Prefer the newest initialization folder if multiple runs exist.
  init_dirs <- list.dirs(tool_dir, full.names = TRUE, recursive = FALSE)
  init_dirs <- init_dirs[grepl("^clone\\d+_rectangle0_w1\\.0$", basename(init_dirs))]
  if (length(init_dirs) == 0) {
    message("  Skip: No clone*_rectangle0_w1.0 folder found.")
    return(NULL)
  }
  init_dir <- init_dirs[order(file.info(init_dirs)$mtime, decreasing = TRUE)][1]
  message("  Using clone folder: ", basename(init_dir))
  
  labels_path <- file.path(init_dir, "clone_labels.tsv")
  gene_path   <- file.path(init_dir, "cnv_genelevel.tsv")
  if (!file.exists(labels_path) || !file.exists(gene_path)) {
    message("  Skip: clone_labels.tsv or cnv_genelevel.tsv missing.")
    return(NULL)
  }
  
  out_log2_file <- file.path(tool_dir, "processed_log2_spots_x_genes.csv")
  out_cnv_file <- file.path(tool_dir, "processed_totalcn_genes_x_spots.csv")
  out_labels_file <- file.path(tool_dir, "processed_clone_labels.csv")
  
  labels <- data.table::fread(labels_path)
  if (ncol(labels) < 2) stop("clone_labels.tsv must have 2 columns: spot, clone_label")
  data.table::setnames(labels, 1:2, c("spot", "clone_label"))
  labels <- labels[, .(
    spot = as.character(spot),
    clone_label = as.character(clone_label)
  )]
  
  gene_dt <- data.table::fread(gene_path)
  if (!"gene" %in% names(gene_dt)) stop("cnv_genelevel.tsv missing 'gene' column")
  
  # Compute total CN per clone: cloneX = cloneX A + cloneX B
  clone_cols <- grep("^clone\\d+\\s+[AB]$", names(gene_dt), value = TRUE)
  clone_ids <- unique(sub("\\s+[AB]$", "", clone_cols))
  if (length(clone_ids) == 0) stop("No clone A/B columns found in cnv_genelevel.tsv")
  clone_ids <- clone_ids[order(as.integer(sub("^clone", "", clone_ids)))]
  
  for (cid in clone_ids) {
    a_col <- paste0(cid, " A")
    b_col <- paste0(cid, " B")
    if (!all(c(a_col, b_col) %in% names(gene_dt))) next
    gene_dt[[cid]] <- gene_dt[[a_col]] + gene_dt[[b_col]]
  }
  
  keep_cols <- c("gene", clone_ids)
  gene_tot <- gene_dt[, ..keep_cols]
  
  labels[, clone_col := paste0("clone", clone_label)]
  missing_cols <- setdiff(unique(labels$clone_col), names(gene_tot))
  if (length(missing_cols) > 0) {
    stop("Missing clone columns in cnv_genelevel.tsv: ", paste(missing_cols, collapse = ", "))
  }
  
  # Broadcast clone CN to spots.
  clone_idx <- match(labels$clone_col, clone_ids)
  clone_mat <- as.matrix(gene_tot[, ..clone_ids])
  storage.mode(clone_mat) <- "numeric"
  gene_x_spot_cn <- clone_mat[, clone_idx, drop = FALSE]
  rownames(gene_x_spot_cn) <- gene_tot$gene
  colnames(gene_x_spot_cn) <- labels$spot
  
  gene_x_spot_dt <- as.data.table(gene_x_spot_cn, keep.rownames = "gene")
  data.table::fwrite(gene_x_spot_dt, out_cnv_file)
  
  labels_out <- labels[, .(spot, clone_label)]
  data.table::fwrite(labels_out, out_labels_file)
  
  # Legacy evaluation file: spots x genes in log2(total_cn / 2).
  mat_log2 <- t(log2((gene_x_spot_cn / 2) + EPS))
  mat_log2_dt <- as.data.table(mat_log2, keep.rownames = "spot")
  data.table::fwrite(mat_log2_dt, out_log2_file)
  
  message("[CalicoST] Done. Saved legacy log2 matrix to: ", out_log2_file)
  message("[CalicoST] Saved total-CN gene x spot matrix to: ", out_cnv_file)
  message("[CalicoST] Saved clone labels to: ", out_labels_file)
}
# ==============================================================================
# Main loop
# ==============================================================================

# Root directory holding the Rep_* batches (override with --root_dir)
MAIN_DIR <- "<DATA_ROOT>/synthetic/processed_data/random_batch_BAF/generate_code/Complete_Batch_v2"
SKIP_COPYVAE <- FALSE
TOOLS_FILTER <- NULL   # NULL = all tools; else character vector of tool names to keep
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if ("--root_dir" %in% args) {
    idx <- which(args == "--root_dir")
    if (length(args) >= idx + 1) {
      MAIN_DIR <- args[idx + 1]
    }
  } else {
    MAIN_DIR <- args[1]
  }
  if ("--skip_copyvae" %in% args) {
    SKIP_COPYVAE <- TRUE
  }
  if ("--tool_out_name" %in% args) {
    idx <- which(args == "--tool_out_name")
    if (length(args) >= idx + 1) {
      TOOL_OUT_NAME <- args[idx + 1]
    }
  }
  if ("--tool_out_suffix" %in% args) {
    idx <- which(args == "--tool_out_suffix")
    if (length(args) >= idx + 1) {
      TOOL_OUT_SUFFIX <- args[idx + 1]
    }
  }
  if ("--infercnv_dir_name" %in% args) {
    idx <- which(args == "--infercnv_dir_name")
    if (length(args) >= idx + 1) {
      INFERCNV_DIR_NAME <- args[idx + 1]
    }
  }
  if ("--tools" %in% args) {
    idx <- which(args == "--tools")
    if (length(args) >= idx + 1) {
      TOOLS_FILTER <- trimws(strsplit(args[idx + 1], ",", fixed = TRUE)[[1]])
      TOOLS_FILTER <- TOOLS_FILTER[nzchar(TOOLS_FILTER)]
      if (length(TOOLS_FILTER) == 0) TOOLS_FILTER <- NULL
    }
  }
}
# Gate: TRUE if a tool should run under the current filter (NULL filter = all).
run_tool <- function(name) is.null(TOOLS_FILTER) || name %in% TOOLS_FILTER
message("[Config] Tools filter = ",
        if (is.null(TOOLS_FILTER)) "<all>" else paste(TOOLS_FILTER, collapse = ", "))
message("[Config] Tool output suffix = ",
        if (nzchar(TOOL_OUT_SUFFIX)) TOOL_OUT_SUFFIX else "<none>")
if (tolower(Sys.getenv("SKIP_COPYVAE", "0")) %in% c("1","true","yes")) {
  SKIP_COPYVAE <- TRUE
}

if (!dir.exists(MAIN_DIR)) stop("MAIN_DIR does not exist!")

# Collect every batch folder starting with "Rep_"
all_batches <- list.dirs(MAIN_DIR, full.names = TRUE, recursive = FALSE)
target_batches <- all_batches[grepl("^Rep_", basename(all_batches))]

message(sprintf("Found %d batches to process.", length(target_batches)))

# Process each batch
for (batch_path in target_batches) {
  message(paste0("\n=== Working on Batch: ", basename(batch_path), " ==="))
  
  # InferCNV
  if (run_tool("InferCNV")) {
    tryCatch({
      run_infercnv_pipeline(batch_path,EPS=EPS)
    }, error = function(e) { message("Error in InferCNV: ", e$message) })
    tryCatch({
      run_infercnv_state_pipeline(batch_path)
    }, error = function(e) { message("Error in InferCNV-state: ", e$message) })
  }
  # copykat
  if (run_tool("CopyKAT")) {
    tryCatch({
      run_copykat_pipeline(batch_path)
    }, error = function(e) { message("Error in Copykat: ", e$message) })
  }
  # infercnvpy
  if (run_tool("InferCNVPy")) {
    tryCatch({
      run_infercnvpy_pipeline(batch_path)
    }, error = function(e) { message("Error in InferCNVpy: ", e$message) })
  }
  # SCEVAN
  if (run_tool("SCEVAN")) {
    tryCatch({
      run_scevan_pipeline(batch_path)
    }, error = function(e) { message("Error in SCEVAN: ", e$message) })
  }
  # CONICSmat
  if (run_tool("CONICSmat")) {
    tryCatch({
      run_conicsmat_pipeline(batch_path)
    }, error = function(e) { message("Error in CONICSmat: ", e$message) })
  }
  # STARCH
  if (run_tool("STARCH")) {
    tryCatch({
      run_starch_pipeline(batch_path, EPS=EPS)
    }, error = function(e) { message("Error in STARCH: ", e$message) })
  }
  # XClone (state-level)
  if (run_tool("XClone")) {
    tryCatch({
      run_xclone_state_pipeline(batch_path)
    }, error = function(e) { message("Error in XClone: ", e$message) })
  }
  # Numbat (state-level)
  if (run_tool("Numbat")) {
    tryCatch({
      run_numbat_state_pipeline(batch_path)
    }, error = function(e) { message("Error in Numbat: ", e$message) })
  }
  # CalicoST (clone -> spot, log2 CN ratio)
  if (run_tool("CalicoST")) {
    tryCatch({
      run_calicost_pipeline(batch_path, EPS=EPS)
    }, error = function(e) { message("Error in CalicoST: ", e$message) })
  }
  
}

message("\nAll Done! Processed files are saved within their respective ", TOOL_OUT_NAME, " folders.")
message("[Config] InferCNV dir name = ", INFERCNV_DIR_NAME)
