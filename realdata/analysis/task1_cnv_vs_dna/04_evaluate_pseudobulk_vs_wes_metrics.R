#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package 'pROC' is required. Please install it before running this script.")
}

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
DEFAULT_OUT_ALL <- NULL  # single-root default resolves to <dataset>/evaluation/pseudobulk_vs_wes_metrics_all_samples.csv
DEFAULT_OUT_ALL_BASENAME <- "pseudobulk_vs_wes_metrics_all_samples.csv"
DEFAULT_OUT_ALL_DATASETS <- normalizePath(
  file.path(DEFAULT_DATASETS_DIR, "..", "results", "task1_cnv_vs_dna", "pseudobulk_vs_wes_metrics_all_datasets.csv"),
  mustWork = FALSE
)
ALLOWED_TOOLS <- c("infercnv", "conicsmat", "copykat", "scevan", "numbat", "xclone","copyvae","infercnvpy","starch","tuscan","calicost")
DEFAULT_PSEUDO_BASENAME <- "pseudobulk_tumor_epi_log2_gene_x_1.csv"
DEFAULT_WES_LOG2_BASENAME <- "wes_gene_log2ratio_matrix.csv"
DEFAULT_WES_CALL_BASENAME <- "wes_gene_call_matrix.csv"
DEFAULT_SAMPLE_OUT_BASENAME <- "pseudobulk_vs_wes_metrics.csv"
DEFAULT_GENE_POS_FILE <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
MIN_POINTS <- 50L

# Map a tool folder to the methods scored from it. infercnv and numbat each yield
# two methods (expr + cnv) from a single folder; every other tool is one method
# keyed on the folder name, reading the default pseudobulk basename.
tool_method_pseudos <- function(tool_name, pseudo_basename) {
  key <- normalize_key(tool_name)
  if (key == "infercnv") {
    return(list(
      list(tool = "infercnv_expr", basename = "pseudobulk_tumor_epi_infercnv_expr_gene_x_1.csv"),
      list(tool = "infercnv_cnv",  basename = "pseudobulk_tumor_epi_infercnv_cnv_gene_x_1.csv")
    ))
  }
  if (key == "numbat") {
    return(list(
      list(tool = "numbat_expr", basename = "pseudobulk_tumor_epi_numbat_expr_gene_x_1.csv"),
      list(tool = "numbat_cnv",  basename = "pseudobulk_tumor_epi_numbat_cnv_gene_x_1.csv")
    ))
  }
  # STARCH is scored on its discrete STATE pseudobulk (mean of -1/0/1 states,
  # centered at 0), the same way infercnv_cnv/numbat_cnv are, rather than on its
  # log2 pseudobulk.
  if (key == "starch") {
    return(list(list(tool = tool_name, basename = "pseudobulk_tumor_epi_state_gene_x_1.csv")))
  }
  list(list(tool = tool_name, basename = pseudo_basename))
}

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out[nzchar(out)]
}

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

choose_col <- function(nms, candidates, required = TRUE, label = "") {
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

normalize_gene_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x
}

normalize_chr <- function(x, want_chr_prefix = TRUE) {
  na_idx <- is.na(x)
  x <- as.character(x)
  x[na_idx] <- NA_character_
  x <- sub("^CHR", "chr", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("M", "MT", "MtDNA", "mitochondria")] <- "M"
  if (want_chr_prefix) ifelse(is.na(x), NA_character_, paste0("chr", x)) else x
}

normalize_call_to_num <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x[x %in% c("GAIN", "AMP", "AMPLIFICATION")] <- "+"
  x[x %in% c("LOSS", "DEL", "DELETION")] <- "-"
  x[x %in% c("NEUTRAL", "DIPLOID")] <- "0"

  out <- rep(NA_integer_, length(x))
  out[x %in% c("-", "-1")] <- -1L
  out[x %in% c("0", "0.0", "+0")] <- 0L
  out[x %in% c("+", "1", "+1")] <- 1L

  suppressWarnings(num <- as.numeric(x))
  out[is.na(out) & !is.na(num) & num < 0] <- -1L
  out[is.na(out) & !is.na(num) & num == 0] <- 0L
  out[is.na(out) & !is.na(num) & num > 0] <- 1L
  out
}

collapse_rows_mean <- function(mat, groups) {
  mat <- as.matrix(mat)
  groups <- as.character(groups)
  sums <- rowsum(mat, group = groups, reorder = FALSE, na.rm = TRUE)
  cnt <- rowsum(!is.na(mat), group = groups, reorder = FALSE)
  sums / pmax(cnt, 1)
}

read_gene_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(fread(path, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L || ncol(df) < 2L) return(NULL)

  genes <- normalize_gene_ids(df[[1]])
  keep <- !is.na(genes) & nzchar(genes)
  if (!any(keep)) return(NULL)

  mat <- as.matrix(df[keep, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- genes[keep]
  colnames(mat) <- colnames(df)[-1]

  if (anyDuplicated(rownames(mat)) > 0L) {
    mat <- collapse_rows_mean(mat, rownames(mat))
  }

  mat
}

build_gene_reference <- function(gene_pos_file) {
  if (!file.exists(gene_pos_file)) stop("gene_pos_file not found: ", gene_pos_file)

  gp <- fread(gene_pos_file, header = FALSE, sep = "\t", quote = "\"")
  if (ncol(gp) < 4L) stop("gene_pos_file must have at least 4 columns: gene, chr, start, end")
  gp <- gp[, 1:4]
  setnames(gp, c("gene", "chr", "start", "end"))

  if (nrow(gp) > 0L && normalize_key(gp$gene[1]) %in% c("gene", "genename", "symbol", "hgncsymbol")) {
    gp <- gp[-1]
  }

  gp[, gene := normalize_gene_ids(gene)]
  gp[, chr := normalize_chr(chr, want_chr_prefix = TRUE)]
  gp[, start := as.integer(round(as.numeric(start)))]
  gp[, end := as.integer(round(as.numeric(end)))]
  gp <- gp[!is.na(gene) & nzchar(gene) & !is.na(chr) & !is.na(start) & !is.na(end) & end >= start]
  gp <- unique(gp, by = "gene")
  setorder(gp, chr, start, end, gene)
  gp[, gene_len := as.numeric(end - start + 1L)]
  gp
}

make_neutral_truth <- function(sample_name, gene_ref) {
  genes <- sort(unique(gene_ref$gene))
  log2_mat <- matrix(0, nrow = length(genes), ncol = 1L, dimnames = list(genes, sample_name))
  call_mat <- matrix(0, nrow = length(genes), ncol = 1L, dimnames = list(genes, sample_name))
  list(log2 = log2_mat, call = call_mat, source = "neutral_no_cnv")
}

project_segments_to_gene_matrices <- function(seg_file, gene_ref) {
  dt <- tryCatch(fread(seg_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)

  col_chr <- choose_col(names(dt), c("CONTIG", "chr", "chromosome", "chrom"), TRUE, "chromosome")
  col_start <- choose_col(names(dt), c("START", "start", "POS_START"), TRUE, "start")
  col_end <- choose_col(names(dt), c("END", "end", "POS_END"), TRUE, "end")
  col_log2 <- choose_col(names(dt), c("MEAN_LOG2_COPY_RATIO", "log2ratio", "log2", "mean_log2_copy_ratio"), TRUE, "log2")
  col_call <- choose_col(names(dt), c("CALL", "call"), TRUE, "call")
  col_wes_sample <- choose_col(names(dt), c("sample", "SAMPLE", "wes_sample", "sample_id", "wes_sample_id"), FALSE, "wes_sample")

  if (is.na(col_wes_sample)) {
    dt[, wes_sample_id := tools::file_path_sans_ext(basename(seg_file))]
  } else {
    dt[, wes_sample_id := as.character(get(col_wes_sample))]
    dt[is.na(wes_sample_id) | wes_sample_id == "", wes_sample_id := tools::file_path_sans_ext(basename(seg_file))]
  }

  dt[, seg_chr := normalize_chr(get(col_chr), want_chr_prefix = TRUE)]
  dt[, seg_start := as.integer(round(as.numeric(get(col_start))))]
  dt[, seg_end := as.integer(round(as.numeric(get(col_end))))]
  dt[, log2ratio := as.numeric(get(col_log2))]
  dt[, call_num := normalize_call_to_num(get(col_call))]
  dt <- dt[!is.na(seg_chr) & !is.na(seg_start) & !is.na(seg_end) & seg_end >= seg_start]
  if (nrow(dt) == 0L) return(NULL)

  seg_iv <- dt[, .(
    sample_id = wes_sample_id,
    chr = seg_chr,
    start = seg_start,
    end = seg_end,
    log2ratio,
    call_num
  )]
  gene_iv <- gene_ref[, .(chr, start, end, gene, gene_len)]
  setkey(seg_iv, chr, start, end)
  setkey(gene_iv, chr, start, end)

  ov <- foverlaps(seg_iv, gene_iv, type = "any", nomatch = NULL)
  if (is.null(ov) || nrow(ov) == 0L) return(NULL)
  ov[, ovlen := pmin(i.end, end) - pmax(i.start, start) + 1L]
  ov <- ov[ovlen > 0L & !is.na(gene_len) & gene_len > 0]
  if (nrow(ov) == 0L) return(NULL)

  agg <- ov[, .(
    log2_weighted = {
      idx <- !is.na(log2ratio)
      if (any(idx)) sum(log2ratio[idx] * ovlen[idx]) / gene_len[1] else NA_real_
    },
    call_weighted = {
      idx <- !is.na(call_num)
      if (any(idx)) sum(call_num[idx] * ovlen[idx]) / gene_len[1] else NA_real_
    }
  ), by = .(sample_id, gene)]

  wes_ids <- sort(unique(agg$sample_id))
  genes <- sort(unique(gene_ref$gene))
  m_log2 <- matrix(NA_real_, nrow = length(genes), ncol = length(wes_ids), dimnames = list(genes, wes_ids))
  m_call <- matrix(NA_real_, nrow = length(genes), ncol = length(wes_ids), dimnames = list(genes, wes_ids))
  r <- match(agg$gene, genes)
  c <- match(agg$sample_id, wes_ids)
  m_log2[cbind(r, c)] <- agg$log2_weighted
  m_call[cbind(r, c)] <- agg$call_weighted

  list(log2 = m_log2, call = m_call, source = seg_file)
}

dataset_dir_from_root <- function(root_dir) {
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  run_root_base <- basename(root_dir)
  if (run_root_base == "run_CNV" && basename(dirname(root_dir)) == "run_links") {
    dirname(dirname(root_dir))
  } else {
    dirname(root_dir)
  }
}

truth_segment_candidates <- function(dataset_dir, dataset_name, sample_name) {
  # Prostate pt10 WGS is low-amplitude: with normal dilution every segment has
  # |log2| < 0.3, so the standard truth is all-neutral (max_F1 = 1/3 for every
  # tool). Score it against the tumor-only renormalized truth instead.
  truth_sub <- if (dataset_name == "metastatic_prostate_ffpe") "wes_cnv_tumor_only" else "wes_cnv"
  truth_dir <- file.path(dataset_dir, "truth_links", truth_sub)
  if (!dir.exists(truth_dir)) return(character(0))

  direct <- file.path(truth_dir, c(
    paste0(sample_name, "_wes_cnv.csv"),
    paste0(sample_name, "_wes_segments_results.csv"),
    sample_name,
    paste0(sample_name, ".csv")
  ))

  mapped <- character(0)
  if (dataset_name == "skin_tumor_ff") {
    patient <- sub("^skin_(P[0-9]+).*$", "\\1", sample_name)
    patient <- sub("^p", "pt", tolower(patient))
    mapped <- file.path(truth_dir, c(
      paste0(patient, "_wes_segments_results.csv"),
      paste0(patient, "_wes_cnv.csv")
    ))
  }

  unique(c(direct, mapped))[file.exists(unique(c(direct, mapped)))]
}

load_truth_matrices <- function(sample_dir, dataset_name, dataset_dir, wes_log2_basename,
                                wes_call_basename, gene_ref) {
  sample_name <- basename(sample_dir)
  wes_log2_file <- file.path(sample_dir, wes_log2_basename)
  wes_call_file <- file.path(sample_dir, wes_call_basename)

  wes_log2 <- read_gene_matrix(wes_log2_file)
  wes_call <- read_gene_matrix(wes_call_file)
  if (!is.null(wes_log2) && nrow(wes_log2) > 0L && ncol(wes_log2) > 0L) {
    return(list(log2 = wes_log2, call = wes_call, source = wes_log2_file))
  }

  if (!(dataset_name %in% c("metastatic_prostate_ffpe", "skin_tumor_ff"))) {
    return(NULL)
  }

  seg_files <- truth_segment_candidates(dataset_dir, dataset_name, sample_name)
  for (seg_file in seg_files) {
    projected <- project_segments_to_gene_matrices(seg_file, gene_ref)
    if (!is.null(projected) && nrow(projected$log2) > 0L && ncol(projected$log2) > 0L) {
      return(projected)
    }
  }

  if (dataset_name == "skin_tumor_ff" && grepl("^skin_P4_", sample_name)) {
    return(make_neutral_truth(sample_name, gene_ref))
  }

  NULL
}

# normalize_call_vector <- function(x) {
#   s_raw <- trimws(as.character(x))
#   s_up <- toupper(s_raw)
# 
#   out <- rep(NA_real_, length(s_raw))
#   out[s_up %in% c("GAIN", "AMP", "AMPLIFICATION", "+", "+1", "1")] <- 1
#   out[s_up %in% c("LOSS", "DEL", "DELETION", "-", "-1")] <- -1
#   out[s_up %in% c("NEUTRAL", "DIPLOID", "0", "+0", "0.0")] <- 0
# 
#   suppressWarnings(num <- as.numeric(s_raw))
#   out[is.na(out) & !is.na(num)] <- num[is.na(out) & !is.na(num)]
# 
#   out
# }

# normalize_call_matrix <- function(mat) {
#   if (is.null(mat)) return(NULL)
#   out <- apply(mat, 2, normalize_call_vector)
#   if (!is.matrix(out)) out <- matrix(out, ncol = 1)
#   rownames(out) <- rownames(mat)
#   colnames(out) <- colnames(mat)
#   out
# }

macro_f1_3class <- function(pred_class, true_class, classes = c(-1L, 0L, 1L)) {
  f1s <- numeric(length(classes))
  names(f1s) <- as.character(classes)

  for (k in seq_along(classes)) {
    cls <- classes[k]
    tp <- sum(pred_class == cls & true_class == cls, na.rm = TRUE)
    fp <- sum(pred_class == cls & true_class != cls, na.rm = TRUE)
    fn <- sum(pred_class != cls & true_class == cls, na.rm = TRUE)

    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1 <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_
    if (is.na(f1)) f1 <- 0
    f1s[k] <- f1
  }

  mean(f1s)
}

get_f1_score_unweighted_log2ratio <- function(pred_score, thr_loss, thr_gain, true_class) {
  pred_class <- ifelse(pred_score < thr_loss, -1L, ifelse(pred_score > thr_gain, 1L, 0L))
  macro_f1_3class(pred_class = pred_class, true_class = true_class)
}

evaluate_optimal_f1_score_log2ratio <- function(pred_score, true_class, baseline = 0, max_candidates = 100) {
  ok <- is.finite(pred_score) & !is.na(true_class)
  pred_score <- pred_score[ok]
  true_class <- as.integer(true_class[ok])

  if (length(pred_score) < MIN_POINTS) {
    return(c(max_f1 = NA_real_, gain_thr = NA_real_, loss_thr = NA_real_))
  }

  scores <- sort(unique(pred_score))

  loss_scores <- scores[scores <= baseline]
  if (length(loss_scores) > max_candidates) {
    loss_scores <- loss_scores[round(seq(1, length(loss_scores), length.out = max_candidates))]
  }
  if (length(loss_scores) == 0L) loss_scores <- baseline
  loss_scores <- unique(c(loss_scores, baseline))

  gain_scores <- scores[scores >= baseline]
  if (length(gain_scores) > max_candidates) {
    gain_scores <- gain_scores[round(seq(1, length(gain_scores), length.out = max_candidates))]
  }
  if (length(gain_scores) == 0L) gain_scores <- baseline
  gain_scores <- unique(c(baseline, gain_scores))

  combis <- expand.grid(loss_thr = loss_scores, gain_thr = gain_scores)
  combis <- combis[combis$loss_thr <= combis$gain_thr, , drop = FALSE]

  combis$f1 <- vapply(seq_len(nrow(combis)), function(i) {
    get_f1_score_unweighted_log2ratio(
      pred_score,
      thr_loss = combis$loss_thr[i],
      thr_gain = combis$gain_thr[i],
      true_class = true_class
    )
  }, numeric(1))

  best <- combis[which.max(combis$f1), , drop = FALSE]
  c(max_f1 = best$f1, gain_thr = best$gain_thr, loss_thr = best$loss_thr)
}

calc_roc_stats <- function(truth_binary, score_vec, score_threshold, direction_gain = TRUE) {
  if (length(unique(truth_binary)) < 2) return(c(auc = NA_real_, pauc = NA_real_))
  if (sum(truth_binary == 1, na.rm = TRUE) == 0) return(c(auc = NA_real_, pauc = NA_real_))

  scores <- if (direction_gain) score_vec else -score_vec
  thr <- if (direction_gain) score_threshold else -score_threshold

  keep <- !is.na(scores) & !is.na(truth_binary)
  scores <- scores[keep]
  y <- truth_binary[keep]

  roc_obj <- suppressMessages(pROC::roc(response = y, predictor = scores, quiet = TRUE))
  auc_val <- as.numeric(pROC::auc(roc_obj))

  pred_pos <- ifelse(scores > thr, 1, 0)
  s_max <- sum(pred_pos == 1 & y == 1) / sum(y == 1)

  pauc_val <- if (!is.na(s_max) && s_max > 0) {
    as.numeric(pROC::auc(
      roc_obj,
      partial.auc = c(0, s_max),
      partial.auc.focus = "sens",
      partial.auc.correct = FALSE
    ))
  } else {
    0
  }

  c(auc = auc_val, pauc = pauc_val)
}

calc_metrics_pseudobulk <- function(pred_score, true_log2, true_call) {
  valid <- is.finite(pred_score) & is.finite(true_log2) & !is.na(true_call)
  p_vec <- pred_score[valid]
  t_log2_vec <- true_log2[valid]
  call_vec <- true_call[valid]

  if (length(p_vec) < MIN_POINTS) {
    return(list(
      Pearson = NA_real_,
      Spearman = NA_real_,
      AUC_gain_trunc = NA_real_,
      AUC_loss_trunc = NA_real_,
      max_F1 = NA_real_,
      valid_points = as.integer(length(p_vec))
    ))
  }

  pear <- if (sd(p_vec) > 0 && sd(t_log2_vec) > 0) cor(p_vec, t_log2_vec,method = 'pearson',use = "complete.obs") else 0
  spear <- if (length(unique(p_vec)) > 1 && length(unique(t_log2_vec)) > 1) {
    cor(p_vec, t_log2_vec, method = "spearman", use = "complete.obs")
  } else {
    0
  }

  gt_gain <- call_vec > 0
  gt_loss <- call_vec < 0
  true_class <- ifelse(call_vec > 0, 1L, ifelse(call_vec < 0, -1L, 0L))

  stats_gain <- calc_roc_stats(as.numeric(gt_gain), p_vec, score_threshold = 0, direction_gain = TRUE)
  stats_loss <- calc_roc_stats(as.numeric(gt_loss), p_vec, score_threshold = 0, direction_gain = FALSE)
  opt <- evaluate_optimal_f1_score_log2ratio(pred_score = p_vec, true_class = true_class, baseline = 0)

  list(
    Pearson = pear,
    Spearman = spear,
    AUC_gain_trunc = unname(stats_gain["pauc"]),
    AUC_loss_trunc = unname(stats_loss["pauc"]),
    max_F1 = unname(opt["max_f1"]),
    valid_points = as.integer(length(p_vec))
  )
}

pick_wes_call_column <- function(wes_id, wes_log2_mat, wes_call_mat) {
  if (is.null(wes_call_mat) || ncol(wes_call_mat) == 0L) return(NA_character_)
  if (wes_id %in% colnames(wes_call_mat)) return(wes_id)

  if (ncol(wes_call_mat) == 1L) return(colnames(wes_call_mat)[1])

  idx <- match(wes_id, colnames(wes_log2_mat))
  if (!is.na(idx) && idx <= ncol(wes_call_mat)) return(colnames(wes_call_mat)[idx])

  NA_character_
}

log_msg <- function(...) cat(sprintf(...), "\n")

process_root <- function(dataset_name, root_dir, out_all_file, pseudo_basename, wes_log2_basename,
                         wes_call_basename, sample_out_basename, samples, strategies, tools,
                         gene_ref) {
  if (!dir.exists(root_dir)) stop("root_dir does not exist: ", root_dir)

  dataset_dir <- dataset_dir_from_root(root_dir)

  # Resolve default out_all_file relative to the dataset folder.
  # Expected layouts: <dataset>/run_links or <dataset>/run_CNV -> merged CSV goes to <dataset>/evaluation/
  if (is.null(out_all_file)) {
    out_all_file <- file.path(dataset_dir, "evaluation", DEFAULT_OUT_ALL_BASENAME)
  }
  dir.create(dirname(out_all_file), showWarnings = FALSE, recursive = TRUE)

  sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[file.info(sample_dirs)$isdir]
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[normalize_key(basename(sample_dirs)) %in% normalize_key(samples)]
  }

  log_msg("[Init] dataset = %s | root_dir = %s", dataset_name, root_dir)
  log_msg("[Init] n_samples = %d", length(sample_dirs))
  log_msg("[Init] pseudo_basename = %s", pseudo_basename)
  log_msg("[Init] wes_log2_basename = %s | wes_call_basename = %s", wes_log2_basename, wes_call_basename)
  log_msg("[Init] dataset_dir = %s", dataset_dir)

  all_rows <- list()
  all_i <- 1L

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    truth <- load_truth_matrices(
      sample_dir = sample_dir,
      dataset_name = dataset_name,
      dataset_dir = dataset_dir,
      wes_log2_basename = wes_log2_basename,
      wes_call_basename = wes_call_basename,
      gene_ref = gene_ref
    )

    if (is.null(truth)) {
      log_msg("[Sample] %s | skip: DNA truth missing/invalid (matrix or truth_links/wes_cnv segment file)", sample_name)
      next
    }
    wes_log2 <- truth$log2
    wes_call <- truth$call

    run_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[file.info(run_dirs)$isdir]
    if (!is.null(strategies)) {
      run_dirs <- run_dirs[normalize_key(basename(run_dirs)) %in% normalize_key(strategies)]
    }

    sample_rows <- list()
    sample_i <- 1L
    log_msg("[Sample] %s | n_wes_samples=%d | n_runs=%d | truth_source=%s",
            sample_name, ncol(wes_log2), length(run_dirs), truth$source)

    for (run_dir in run_dirs) {
      run_name <- basename(run_dir)
      tool_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
      tool_dirs <- tool_dirs[file.info(tool_dirs)$isdir]
      # Whitelist: only evaluate tools in ALLOWED_TOOLS (case-insensitive match on folder name).
      # This also naturally excludes figures/tables/*-old and any other non-tool subdirs.
      tool_dirs <- tool_dirs[normalize_key(basename(tool_dirs)) %in% normalize_key(ALLOWED_TOOLS)]
      if (!is.null(tools)) {
        tool_dirs <- tool_dirs[normalize_key(basename(tool_dirs)) %in% normalize_key(tools)]
      }

      for (tool_dir in tool_dirs) {
        folder_name <- basename(tool_dir)
        for (method in tool_method_pseudos(folder_name, pseudo_basename)) {
        tool_name <- method$tool
        pseudo_file <- file.path(tool_dir, method$basename)

        # Numbat writes no numbat_cnv pseudobulk when no CNV survives its LLR
        # filter. That is a valid flat (all-neutral) call, so score it as
        # Spearman 0 instead of dropping it as pseudo_missing.
        numbat_log <- file.path(tool_dir, "log.txt")
        if (!file.exists(pseudo_file) && tool_name == "numbat_cnv" && file.exists(numbat_log) &&
            any(grepl("No CNV remains after filtering by LLR", readLines(numbat_log, warn = FALSE), fixed = TRUE))) {
          sample_rows[[sample_i]] <- data.table(
            sample = sample_name,
            run = run_name,
            tool = tool_name,
            wes_sample_id = NA_character_,
            n_genes_overlap = 0L,
            valid_points = 0L,
            Pearson = NA_real_,
            Spearman = 0,
            AUC_gain_trunc = NA_real_,
            AUC_loss_trunc = NA_real_,
            max_F1 = NA_real_,
            status = "ok",
            pseudo_file = pseudo_file,
            call_source = "numbat_no_cnv_detected"
          )
          sample_i <- sample_i + 1L
          next
        }

        if (!file.exists(pseudo_file)) {
          sample_rows[[sample_i]] <- data.table(
            sample = sample_name,
            run = run_name,
            tool = tool_name,
            wes_sample_id = NA_character_,
            n_genes_overlap = 0L,
            valid_points = 0L,
            Pearson = NA_real_,
            Spearman = NA_real_,
            AUC_gain_trunc = NA_real_,
            AUC_loss_trunc = NA_real_,
            max_F1 = NA_real_,
            status = "pseudo_missing",
            pseudo_file = pseudo_file,
            call_source = NA_character_
          )
          sample_i <- sample_i + 1L
          next
        }

        pred_mat <- read_gene_matrix(pseudo_file)
        if (is.null(pred_mat) || nrow(pred_mat) == 0L || ncol(pred_mat) == 0L) {
          sample_rows[[sample_i]] <- data.table(
            sample = sample_name,
            run = run_name,
            tool = tool_name,
            wes_sample_id = NA_character_,
            n_genes_overlap = 0L,
            valid_points = 0L,
            Pearson = NA_real_,
            Spearman = NA_real_,
            AUC_gain_trunc = NA_real_,
            AUC_loss_trunc = NA_real_,
            max_F1 = NA_real_,
            status = "pseudo_invalid",
            pseudo_file = pseudo_file,
            call_source = NA_character_
          )
          sample_i <- sample_i + 1L
          next
        }

        pred_vec <- pred_mat[, 1]
        names(pred_vec) <- rownames(pred_mat)

        for (wes_id in colnames(wes_log2)) {
          gt_log2_vec <- wes_log2[, wes_id]
          names(gt_log2_vec) <- rownames(wes_log2)

          call_col <- pick_wes_call_column(wes_id, wes_log2, wes_call)
          if (!is.na(call_col)) {
            gt_call_vec <- wes_call[, call_col]
            names(gt_call_vec) <- rownames(wes_call)
            call_source <- "wes_call_matrix"
          } else {
            gt_call_vec <- ifelse(gt_log2_vec > 0, 1, ifelse(gt_log2_vec < 0, -1, 0))
            call_source <- "fallback_sign_wes_log2"
          }

          common_genes <- intersect(names(pred_vec), intersect(names(gt_log2_vec), names(gt_call_vec)))
          n_overlap <- length(common_genes)

          if (n_overlap == 0L) {
            sample_rows[[sample_i]] <- data.table(
              sample = sample_name,
              run = run_name,
              tool = tool_name,
              wes_sample_id = wes_id,
              n_genes_overlap = 0L,
              valid_points = 0L,
              Pearson = NA_real_,
              Spearman = NA_real_,
              AUC_gain_trunc = NA_real_,
              AUC_loss_trunc = NA_real_,
              max_F1 = NA_real_,
              status = "no_gene_overlap",
              pseudo_file = pseudo_file,
              call_source = call_source
            )
            sample_i <- sample_i + 1L
            next
          }

          p <- as.numeric(pred_vec[common_genes])
          t_log2 <- as.numeric(gt_log2_vec[common_genes])
          t_call <- as.numeric(gt_call_vec[common_genes])
          mt <- calc_metrics_pseudobulk(p, t_log2, t_call)

          status <- if (mt$valid_points < MIN_POINTS) "too_few_valid_points" else "ok"

          sample_rows[[sample_i]] <- data.table(
            sample = sample_name,
            run = run_name,
            tool = tool_name,
            wes_sample_id = wes_id,
            n_genes_overlap = as.integer(n_overlap),
            valid_points = as.integer(mt$valid_points),
            Pearson = as.numeric(mt$Pearson),
            Spearman = as.numeric(mt$Spearman),
            AUC_gain_trunc = as.numeric(mt$AUC_gain_trunc),
            AUC_loss_trunc = as.numeric(mt$AUC_loss_trunc),
            max_F1 = as.numeric(mt$max_F1),
            status = status,
            pseudo_file = pseudo_file,
            call_source = call_source
          )
          sample_i <- sample_i + 1L
        }
        }
      }
    }

    if (length(sample_rows) > 0L) {
      sample_dt <- rbindlist(sample_rows, fill = TRUE)
      sample_dt[, dataset := dataset_name]
      setcolorder(sample_dt, c("dataset", setdiff(names(sample_dt), "dataset")))
      out_sample <- file.path(sample_dir, sample_out_basename)
      fwrite(sample_dt, out_sample)
      log_msg("[Save] %s | rows=%d | %s", sample_name, nrow(sample_dt), out_sample)

      all_rows[[all_i]] <- sample_dt
      all_i <- all_i + 1L
    }
  }

  if (length(all_rows) == 0L) {
    log_msg("[Done] %s | No evaluation rows generated.", dataset_name)
    return(NULL)
  }

  all_dt <- rbindlist(all_rows, fill = TRUE)
  fwrite(all_dt, out_all_file)

  log_msg("[Done] %s | total_rows=%d | ok=%d | too_few=%d | pseudo_missing=%d | no_overlap=%d",
          dataset_name,
          nrow(all_dt),
          sum(all_dt$status == "ok"),
          sum(all_dt$status == "too_few_valid_points"),
          sum(all_dt$status == "pseudo_missing"),
          sum(all_dt$status == "no_gene_overlap"))
  log_msg("[Done] all-sample result: %s", out_all_file)

  all_dt
}

main <- function() {
  root_dir <- DEFAULT_ROOT_DIR
  datasets_dir <- DEFAULT_DATASETS_DIR
  out_all_file <- DEFAULT_OUT_ALL
  pseudo_basename <- DEFAULT_PSEUDO_BASENAME
  wes_log2_basename <- DEFAULT_WES_LOG2_BASENAME
  wes_call_basename <- DEFAULT_WES_CALL_BASENAME
  sample_out_basename <- DEFAULT_SAMPLE_OUT_BASENAME
  gene_pos_file <- DEFAULT_GENE_POS_FILE
  datasets <- NULL
  samples <- NULL
  strategies <- NULL
  tools <- NULL

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
      if (key == "--out_all_file" && i + 1L <= length(args)) {
        out_all_file <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--pseudo_basename" && i + 1L <= length(args)) {
        pseudo_basename <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--wes_log2_basename" && i + 1L <= length(args)) {
        wes_log2_basename <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--wes_call_basename" && i + 1L <= length(args)) {
        wes_call_basename <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--sample_out_basename" && i + 1L <= length(args)) {
        sample_out_basename <- args[i + 1L]; i <- i + 2L; next
      }
      if (key == "--gene_pos_file" && i + 1L <= length(args)) {
        gene_pos_file <- args[i + 1L]; i <- i + 2L; next
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
      i <- i + 1L
    }
  }

  roots <- if (!is.null(root_dir) && nzchar(root_dir)) {
    data.table(dataset = dataset_name_from_root(root_dir), root_dir = normalizePath(root_dir, mustWork = TRUE))
  } else {
    discover_dataset_roots(datasets_dir = datasets_dir, datasets = datasets)
  }
  if (nrow(roots) == 0L) {
    log_msg("[Done] No dataset roots found.")
    return(invisible(NULL))
  }

  gene_ref <- build_gene_reference(gene_pos_file)
  log_msg("[Init] gene_pos_file = %s | genes = %d", gene_pos_file, nrow(gene_ref))

  all_rows <- list()
  for (i in seq_len(nrow(roots))) {
    dataset_out <- if (nrow(roots) == 1L) out_all_file else NULL
    res <- process_root(
      dataset_name = roots$dataset[i],
      root_dir = roots$root_dir[i],
      out_all_file = dataset_out,
      pseudo_basename = pseudo_basename,
      wes_log2_basename = wes_log2_basename,
      wes_call_basename = wes_call_basename,
      sample_out_basename = sample_out_basename,
      samples = samples,
      strategies = strategies,
      tools = tools,
      gene_ref = gene_ref
    )
    if (!is.null(res)) all_rows[[length(all_rows) + 1L]] <- res
  }

  if (length(all_rows) > 0L && nrow(roots) > 1L) {
    all_dt <- rbindlist(all_rows, fill = TRUE)
    combined_out <- if (!is.null(out_all_file)) out_all_file else DEFAULT_OUT_ALL_DATASETS
    dir.create(dirname(combined_out), showWarnings = FALSE, recursive = TRUE)
    fwrite(all_dt, combined_out)
    log_msg("[Done] all_datasets=%d | total_rows=%d | ok=%d | combined_result=%s",
            length(unique(all_dt$dataset)),
            nrow(all_dt),
            sum(all_dt$status == "ok"),
            combined_out)
  }
}

main()
