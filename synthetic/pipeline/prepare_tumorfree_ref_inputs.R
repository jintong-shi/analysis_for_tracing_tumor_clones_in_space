#!/usr/bin/env Rscript
# Standardized inputs for the reference-dependent tools, taking the normal
# reference from a separate tumor-free slide (--ref_dir) instead of the same slide.
# Writes two reference sets per replicate:
#   Standardized_Inputs_TumorFreeNormal   pure-normal spots (Tumor_prop < cutoff)
#   Standardized_Inputs_TumorFreeTME      TME-matched spots
suppressPackageStartupMessages({
  library(arrow)
  library(Matrix)
})

parse_args <- function(args) {
  opt <- list(
    rep_root = NULL,
    ref_dir = NULL,
    batch_root = NULL,
    reps = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20",
    tumor_prop_cutoff = "0.2",
    normal_std_in = "Standardized_Inputs_TumorFreeNormal",
    tme_std_in = "Standardized_Inputs_TumorFreeTME",
    force = "FALSE"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- substring(key, 3)
    if (!name %in% names(opt)) stop("Unknown argument: ", key)
    if (i == length(args)) stop("Missing value after ", key)
    opt[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  missing <- names(opt)[vapply(opt, is.null, logical(1))]
  if (length(missing) > 0) stop("Missing required args: ", paste(missing, collapse = ", "))
  opt$reps <- as.integer(strsplit(opt$reps, ",", fixed = TRUE)[[1]])
  opt$tumor_prop_cutoff <- as.numeric(opt$tumor_prop_cutoff)
  opt$force <- toupper(opt$force) %in% c("TRUE", "T", "1", "YES")
  opt
}

resolve_rep_dir <- function(root, rep_id) {
  candidates <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  candidates <- candidates[grepl(paste0("^Rep_", rep_id, "_"), basename(candidates))]
  if (length(candidates) != 1) {
    stop("Expected exactly one directory for Rep_", rep_id, ", found ", length(candidates))
  }
  candidates[[1]]
}

read_spatial_counts <- function(expr_file) {
  if (!file.exists(expr_file)) stop("Missing expression parquet: ", expr_file)
  df <- arrow::read_parquet(expr_file)
  if ("index" %in% colnames(df)) {
    genes <- as.character(df$index)
    df$index <- NULL
  } else {
    genes <- as.character(df[[1]])
    df[[1]] <- NULL
  }
  mat <- as.matrix(df)
  rownames(mat) <- genes
  Matrix(mat, sparse = TRUE)
}

read_spatial_table <- function(path) {
  if (!file.exists(path)) stop("Missing parquet: ", path)
  df <- as.data.frame(arrow::read_parquet(path), check.names = FALSE)
  colnames(df)[1] <- "cell_id"
  df$cell_id <- as.character(df$cell_id)
  df
}

read_sc_meta <- function(path) {
  if (!file.exists(path)) stop("Missing single-cell metadata: ", path)
  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1] <- "cell_id"
  df$cell_id <- as.character(df$cell_id)
  df$cell_type <- as.character(df$cell_type)
  df
}

parse_cell_by_spot <- function(path, sc_ids) {
  if (!file.exists(path)) stop("Missing cell-by-spot map: ", path)
  lines <- readLines(path, warn = FALSE)
  out <- vector("list", length(lines))
  sc_prefixed <- any(grepl("^Cell_", head(sc_ids, min(5, length(sc_ids)))))
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[[i]], "\\s*:\\s*")[[1]]
    if (length(parts) < 2) next
    spot <- trimws(parts[[1]])
    cells <- trimws(unlist(strsplit(parts[[2]], ",", fixed = TRUE)))
    cells <- cells[nzchar(cells)]
    if (length(cells) > 0 && sc_prefixed && !any(grepl("^Cell_", head(cells, min(5, length(cells)))))) {
      cells <- paste0("Cell_", cells)
    }
    out[[i]] <- list(spot = spot, cells = cells)
  }
  Filter(Negate(is.null), out)
}

spot_celltype_props <- function(run_dir) {
  annot <- read_spatial_table(file.path(run_dir, "synthetic", "spatial_data", "st-annotation.parquet"))
  sc_meta <- read_sc_meta(file.path(run_dir, "synthetic", "single_cell_data", "sc-meta_data.tsv"))
  parsed <- parse_cell_by_spot(
    file.path(run_dir, "synthetic", "spatial_data", "cell_by_spot.txt"),
    sc_meta$cell_id
  )
  type_by_cell <- setNames(sc_meta$cell_type, sc_meta$cell_id)
  all_types <- sort(unique(sc_meta$cell_type))
  prop_mat <- matrix(0, nrow = length(parsed), ncol = length(all_types),
                     dimnames = list(vapply(parsed, `[[`, character(1), "spot"), all_types))
  for (i in seq_along(parsed)) {
    types <- type_by_cell[parsed[[i]]$cells]
    types <- types[!is.na(types)]
    if (length(types) == 0) next
    tab <- table(types) / length(types)
    prop_mat[i, names(tab)] <- as.numeric(tab)
  }
  props <- data.frame(cell_id = rownames(prop_mat), prop_mat, row.names = NULL, check.names = FALSE)
  tumor_cols <- grep("Tumor", colnames(props), ignore.case = TRUE, value = TRUE)
  props$Tumor_prop <- if (length(tumor_cols) > 0) rowSums(props[, tumor_cols, drop = FALSE]) else 0
  props$Hepatocyte_prop <- if ("Hepatocyte" %in% colnames(props)) props$Hepatocyte else 0
  tme_cols <- intersect(c("Fibroblast", "Endothelial", "T", "T_reg", "B", "NK", "Macrophage", "Plasma"), colnames(props))
  props$TME_score <- if (length(tme_cols) > 0) rowSums(props[, tme_cols, drop = FALSE]) else 0
  merge(annot, props, by = "cell_id", all.x = TRUE, sort = FALSE)
}

read_baf <- function(baf_dir) {
  barcodes <- readLines(file.path(baf_dir, "barcodes.tsv"), warn = FALSE)
  ad <- readMM(file.path(baf_dir, "AD.mtx"))
  dp <- readMM(file.path(baf_dir, "DP.mtx"))
  rownames(ad) <- barcodes
  rownames(dp) <- barcodes
  list(
    ad = as(ad, "dgCMatrix"),
    dp = as(dp, "dgCMatrix"),
    barcodes = barcodes,
    snps_file = file.path(baf_dir, "snps.tsv"),
    snp_info_file = file.path(baf_dir, "snp_info.csv")
  )
}

link_or_copy <- function(src, dst) {
  if (file.exists(dst) || file.symlink(src, dst)) return(invisible(TRUE))
  file.copy(src, dst, overwrite = TRUE)
}

write_expr_parquet <- function(counts, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  dense <- as.matrix(counts)
  df <- data.frame(index = rownames(dense), dense, check.names = FALSE)
  arrow::write_parquet(df, out_file)
  invisible(TRUE)
}

write_meta_parquet <- function(meta_df, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(meta_df, out_file)
}

write_numbat_alleles <- function(ad, dp, snp_info, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  dp_s <- summary(dp)
  if (nrow(dp_s) == 0) stop("No nonzero DP entries available for Numbat allele file")
  ad_vals <- as.numeric(ad[cbind(dp_s$i, dp_s$j)])
  parts <- strsplit(as.character(snp_info$snp_id), "_", fixed = TRUE)
  ref <- vapply(parts, function(x) if (length(x) >= 3) x[[3]] else "N", character(1))
  alt <- vapply(parts, function(x) if (length(x) >= 4) x[[4]] else "N", character(1))
  snp_meta <- data.frame(
    snp_id = as.character(snp_info$snp_id),
    CHROM = as.character(snp_info$chrom),
    POS = as.integer(snp_info$pos),
    cM = as.numeric(snp_info$pos) / 1e6,
    REF = ref,
    ALT = alt,
    gene = as.character(snp_info$gene_id),
    stringsAsFactors = FALSE
  )
  con <- gzfile(out_file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("cell\tsnp_id\tCHROM\tPOS\tcM\tREF\tALT\tAD\tDP\tGT\tgene", con)
  chunk <- 1000000L
  starts <- seq.int(1L, nrow(dp_s), by = chunk)
  for (start in starts) {
    end <- min(start + chunk - 1L, nrow(dp_s))
    idx <- start:end
    m <- snp_meta[dp_s$j[idx], , drop = FALSE]
    out <- data.frame(
      cell = rownames(dp)[dp_s$i[idx]],
      snp_id = m$snp_id,
      CHROM = m$CHROM,
      POS = m$POS,
      cM = m$cM,
      REF = m$REF,
      ALT = m$ALT,
      AD = as.integer(round(ad_vals[idx])),
      DP = as.integer(round(dp_s$x[idx])),
      GT = "0|1",
      gene = m$gene,
      stringsAsFactors = FALSE
    )
    write.table(out, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

write_strategy_inputs <- function(rep_dir, shadow_rep_dir, ref_dir, strategy, std_name,
                                  ref_meta, ref_mat, ref_baf, tumor_prop_cutoff, force) {
  std_in <- file.path(shadow_rep_dir, std_name)
  done_file <- file.path(std_in, "prep_done.txt")
  if (file.exists(done_file) && !force) {
    message("[skip] ", basename(rep_dir), " ", strategy, " inputs already prepared")
    return(data.frame(rep = basename(rep_dir), strategy = strategy, std_in = std_in, status = "skipped"))
  }
  if (dir.exists(std_in) && force) unlink(std_in, recursive = TRUE)
  dir.create(std_in, recursive = TRUE, showWarnings = FALSE)

  target_meta <- spot_celltype_props(rep_dir)
  target_meta$Tumor_prop[is.na(target_meta$Tumor_prop)] <- 0
  tumor_spots <- target_meta$cell_id[target_meta$Tumor_prop > tumor_prop_cutoff]
  if (length(tumor_spots) == 0) stop("No tumor spots selected for ", basename(rep_dir))

  if (strategy == "normal") {
    ref_spots <- ref_meta$cell_id[ref_meta$Tumor_prop == 0 & ref_meta$Hepatocyte_prop > 0.5]
  } else if (strategy == "tme") {
    ref_spots <- ref_meta$cell_id[ref_meta$Tumor_prop == 0 & ref_meta$TME_score > 0.5]
  } else {
    stop("Unknown strategy: ", strategy)
  }
  if (length(ref_spots) == 0) stop("No tumor-free reference spots selected for strategy: ", strategy)

  target_mat <- read_spatial_counts(file.path(rep_dir, "synthetic", "spatial_data", "st-expression.parquet"))
  common_genes <- rownames(target_mat)[rownames(target_mat) %in% rownames(ref_mat)]
  target_counts <- target_mat[common_genes, tumor_spots, drop = FALSE]
  ref_counts <- ref_mat[common_genes, ref_spots, drop = FALSE]
  ref_barcodes <- paste0("TumorFree__", ref_spots)
  colnames(ref_counts) <- ref_barcodes
  counts <- cbind(target_counts, ref_counts)

  ann <- data.frame(
    cell_id = c(tumor_spots, ref_barcodes),
    group = c(rep("Tumor", length(tumor_spots)), rep("Normal", length(ref_barcodes))),
    stringsAsFactors = FALSE
  )
  refs <- ref_barcodes

  dir_infer <- file.path(std_in, "InferCNV")
  dir.create(dir_infer, recursive = TRUE, showWarnings = FALSE)
  saveRDS(counts, file.path(dir_infer, "counts.rds"))
  write.table(setNames(data.frame(ann$group, row.names = ann$cell_id), "group"),
              file.path(dir_infer, "annotations.txt"), sep = "\t", quote = FALSE)

  for (tool in c("CopyKAT", "SCEVAN", "CONICSmat")) {
    tool_dir <- file.path(std_in, tool)
    dir.create(tool_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(dir_infer, "counts.rds"), file.path(tool_dir, "counts.rds"), overwrite = TRUE)
    saveRDS(refs, file.path(tool_dir, "references.rds"))
  }

  nb_dir <- file.path(std_in, "Numbat")
  dir.create(nb_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dir_infer, "counts.rds"), file.path(nb_dir, "counts.rds"), overwrite = TRUE)
  write.table(ann, file.path(nb_dir, "annotations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  ipy_dir <- file.path(std_in, "InferCNVPy")
  dir.create(ipy_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(refs, file.path(ipy_dir, "ref_spots.txt"))

  cs_dir <- file.path(std_in, "ClonalScope")
  dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
  writeMM(counts, file.path(cs_dir, "counts.mtx"))
  writeLines(rownames(counts), file.path(cs_dir, "genes.tsv"))
  writeLines(colnames(counts), file.path(cs_dir, "barcodes.tsv"))
  cs_ann <- data.frame(group = tolower(ann$group), row.names = ann$cell_id)
  write.table(cs_ann, file.path(cs_dir, "annotations.txt"), sep = "\t", quote = FALSE)

  target_baf <- read_baf(file.path(rep_dir, "baf"))
  ad_target <- target_baf$ad[tumor_spots, , drop = FALSE]
  dp_target <- target_baf$dp[tumor_spots, , drop = FALSE]
  ad_ref <- ref_baf$ad[ref_spots, , drop = FALSE]
  dp_ref <- ref_baf$dp[ref_spots, , drop = FALSE]
  rownames(ad_ref) <- ref_barcodes
  rownames(dp_ref) <- ref_barcodes
  ad_combined <- rbind(ad_target, ad_ref)
  dp_combined <- rbind(dp_target, dp_ref)

  xclone_dir <- file.path(std_in, "XClone")
  xclone_baf_dir <- file.path(xclone_dir, "BAF")
  xclone_anno_dir <- file.path(xclone_dir, "Annotation")
  xclone_rdr_dir <- file.path(xclone_dir, "RDR")
  dir.create(xclone_baf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(xclone_anno_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(xclone_rdr_dir, recursive = TRUE, showWarnings = FALSE)
  writeMM(ad_combined, file.path(xclone_baf_dir, "AD.mtx"))
  writeMM(dp_combined, file.path(xclone_baf_dir, "DP.mtx"))
  writeLines(rownames(ad_combined), file.path(xclone_baf_dir, "barcodes.tsv"))
  file.copy(target_baf$snps_file, file.path(xclone_baf_dir, "snps.tsv"), overwrite = TRUE)
  file.copy(target_baf$snp_info_file, file.path(xclone_baf_dir, "snp_info.csv"), overwrite = TRUE)
  write.csv(ann, file.path(xclone_dir, "annotation.csv"), row.names = FALSE, quote = FALSE)
  write.csv(ann, file.path(xclone_anno_dir, "annotation.csv"), row.names = FALSE, quote = FALSE)

  snp_info <- read.csv(target_baf$snp_info_file, stringsAsFactors = FALSE, check.names = FALSE)
  write_numbat_alleles(ad_combined, dp_combined, snp_info, file.path(nb_dir, "alleles.tsv.gz"))

  combined_dir <- file.path(std_in, "Combined_Data", "spatial_data")
  expr_parquet <- file.path(combined_dir, "st-expression.parquet")
  annot_parquet <- file.path(combined_dir, "st-annotation.parquet")
  meta_parquet <- file.path(combined_dir, "st-meta.parquet")
  write_expr_parquet(counts, expr_parquet)
  write_meta_parquet(data.frame(index = ann$cell_id, annotation = ann$group, check.names = FALSE), annot_parquet)

  target_spatial_meta <- read_spatial_table(file.path(rep_dir, "synthetic", "spatial_data", "st-meta.parquet"))
  ref_spatial_meta <- read_spatial_table(file.path(ref_dir, "synthetic", "spatial_data", "st-meta.parquet"))
  target_spatial_meta <- target_spatial_meta[match(tumor_spots, target_spatial_meta$cell_id), , drop = FALSE]
  ref_spatial_meta <- ref_spatial_meta[match(ref_spots, ref_spatial_meta$cell_id), , drop = FALSE]
  ref_spatial_meta$cell_id <- ref_barcodes
  combined_meta <- rbind(target_spatial_meta, ref_spatial_meta)
  colnames(combined_meta)[1] <- "index"
  write_meta_parquet(combined_meta, meta_parquet)

  dir.create(file.path(ipy_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  link_or_copy(combined_dir, file.path(ipy_dir, "data", "spatial_data"))
  link_or_copy(expr_parquet, file.path(xclone_rdr_dir, "st-expression.parquet"))
  link_or_copy(meta_parquet, file.path(xclone_rdr_dir, "st-meta.parquet"))

  target_meta$group <- ifelse(target_meta$Tumor_prop > tumor_prop_cutoff, "Tumor", "Ignore")
  ref_meta_out <- ref_meta[match(ref_spots, ref_meta$cell_id), , drop = FALSE]
  ref_meta_out$cell_id <- ref_barcodes
  ref_meta_out$group <- "Normal"
  metadata_out <- rbind.fill.simple(target_meta, ref_meta_out)
  write.csv(metadata_out, file.path(std_in, "Metadata_cellprop.csv"), row.names = FALSE)
  writeLines(refs, file.path(std_in, paste0(strategy, "_ref_spots.txt")))

  writeLines(c(
    paste0("rep=", basename(rep_dir)),
    paste0("strategy=", strategy),
    paste0("tumor_spots=", length(tumor_spots)),
    paste0("ref_spots=", length(ref_spots)),
    paste0("genes=", nrow(counts))
  ), done_file)

  data.frame(
    rep = basename(rep_dir),
    strategy = strategy,
    std_in = std_in,
    n_genes = nrow(counts),
    n_tumor_spots = length(tumor_spots),
    n_ref_spots = length(ref_spots),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

rbind.fill.simple <- function(x, y) {
  all_cols <- union(colnames(x), colnames(y))
  for (nm in setdiff(all_cols, colnames(x))) x[[nm]] <- NA
  for (nm in setdiff(all_cols, colnames(y))) y[[nm]] <- NA
  rbind(x[, all_cols, drop = FALSE], y[, all_cols, drop = FALSE])
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$batch_root, recursive = TRUE, showWarnings = FALSE)

  message("Loading tumor-free reference once: ", opt$ref_dir)
  ref_meta <- spot_celltype_props(opt$ref_dir)
  ref_meta$Tumor_prop[is.na(ref_meta$Tumor_prop)] <- 0
  ref_meta$Hepatocyte_prop[is.na(ref_meta$Hepatocyte_prop)] <- 0
  ref_meta$TME_score[is.na(ref_meta$TME_score)] <- 0
  ref_mat <- read_spatial_counts(file.path(opt$ref_dir, "synthetic", "spatial_data", "st-expression.parquet"))
  ref_baf <- read_baf(file.path(opt$ref_dir, "baf"))

  summaries <- list()
  for (rep_id in opt$reps) {
    rep_dir <- resolve_rep_dir(opt$rep_root, rep_id)
    shadow_rep_dir <- file.path(opt$batch_root, basename(rep_dir))
    if (!dir.exists(shadow_rep_dir)) stop("Missing shadow replicate directory: ", shadow_rep_dir)
    message("\n=== ", basename(rep_dir), " ===")
    summaries[[length(summaries) + 1L]] <- write_strategy_inputs(
      rep_dir, shadow_rep_dir, opt$ref_dir, "normal", opt$normal_std_in,
      ref_meta, ref_mat, ref_baf, opt$tumor_prop_cutoff, opt$force
    )
    summaries[[length(summaries) + 1L]] <- write_strategy_inputs(
      rep_dir, shadow_rep_dir, opt$ref_dir, "tme", opt$tme_std_in,
      ref_meta, ref_mat, ref_baf, opt$tumor_prop_cutoff, opt$force
    )
    gc()
  }
  summary_df <- do.call(rbind, summaries)
  write.csv(summary_df, file.path(opt$batch_root, "tumorfree_ref_input_summary.csv"), row.names = FALSE)
  print(summary_df)
}

main()
