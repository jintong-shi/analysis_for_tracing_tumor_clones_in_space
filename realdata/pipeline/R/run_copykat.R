#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(copykat)
  library(readr)
})

option_list <- list(
  make_option("--visium_folder",      type="character"),
  make_option("--output_dir",         type="character"),
  make_option("--keep_barcodes_file", type="character"),
  make_option("--ref_barcodes_file",  type="character", default=NULL),
  make_option("--threads",            type="integer",   default=4)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("visium_folder","output_dir","keep_barcodes_file"))
  if (is.null(opt[[k]])) stop("--", k, " is required")

## -- 统一创建并切到输出目录（绝对路径） --
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
outdir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
message("[copykat] outdir=", outdir)
if (!dir.exists(outdir)) stop("Output dir not found after dir.create: ", outdir)
old_wd <- getwd()
setwd(outdir)
on.exit(setwd(old_wd), add = TRUE)

## -- 读 10x 计数 --
get_counts_from_visium <- function(visium_folder) {
  mtx_dirs <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix"),
    file.path(visium_folder, "filtered_feature_bc_matrix"),
    visium_folder
  )
  has_mtx <- function(d) file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))
  mtx_dir <- mtx_dirs[file.exists(mtx_dirs) & vapply(mtx_dirs, has_mtx, logical(1))]
  if (length(mtx_dir) > 0) { message("[copykat] Using 10x MTX: ", mtx_dir[1]); return(Seurat::Read10X(mtx_dir[1])) }
  h5_files <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_files <- h5_files[file.exists(h5_files)]
  if (length(h5_files) > 0) { message("[copykat] Using 10x H5: ", h5_files[1]); return(Seurat::Read10X_h5(h5_files[1])) }
  stop("Neither 10x MEX dir nor filtered_feature_bc_matrix.h5 found: ", visium_folder)
}

trim_nonempty <- function(x) unique(x[nchar(x) > 0])

counts <- get_counts_from_visium(opt$visium_folder)

## -- keep/ref 匹配 & 审计 --
keep_raw <- trim_nonempty(readr::read_lines(opt$keep_barcodes_file))
ref_raw  <- character(0)
if (!is.null(opt$ref_barcodes_file) && nzchar(opt$ref_barcodes_file) && file.exists(opt$ref_barcodes_file)) {
  ref_raw <- trim_nonempty(readr::read_lines(opt$ref_barcodes_file))
}
keep_matched <- intersect(keep_raw, colnames(counts))
ref_matched  <- intersect(ref_raw,  colnames(counts))

writeLines(c(
  sprintf("keep_path: %s", normalizePath(opt$keep_barcodes_file, winslash="/", mustWork=FALSE)),
  sprintf(" ref_path: %s", if (is.null(opt$ref_barcodes_file)) "(none)" else normalizePath(opt$ref_barcodes_file, winslash="/", mustWork=FALSE)),
  sprintf("keep_raw: %d", length(keep_raw)),
  sprintf(" ref_raw: %d", length(ref_raw)),
  sprintf("keep_matched: %d", length(keep_matched)),
  sprintf(" ref_matched: %d", length(ref_matched))
), "REF_AUDIT.txt")
writeLines(keep_matched, "keep_used.txt")
writeLines(ref_matched,  "ref_used.txt")
message(sprintf("[copykat][audit] keep_raw=%d keep_matched=%d; ref_raw=%d ref_matched=%d",
                length(keep_raw), length(keep_matched), length(ref_raw), length(ref_matched)))

## -- 选择列 = union(ref, keep)，保持原始顺序 --
sel <- union(ref_matched, keep_matched)
sel <- colnames(counts)[colnames(counts) %in% sel]
if (length(sel) == 0) stop("No cells selected after union(ref, keep).")
expr <- counts[, sel, drop = FALSE]

## -- 跑 copykat（有 ref 就传 norm.cell.names） --
set.seed(1234)
fit <- try(
  if (length(ref_matched) > 0) {
    copykat::copykat(rawmat = expr, norm.cell.names = ref_matched, n.cores = max(1, opt$threads))
  } else {
    copykat::copykat(rawmat = expr, n.cores = max(1, opt$threads))
  },
  silent = TRUE
)

## -- 保存结果：在 outdir 下用相对路径，并做诊断 --
tryCatch({
  saveRDS(fit, file = "copykat_results.rds")
}, error = function(e) {
  message("[copykat] saveRDS failed: ", conditionMessage(e))
  message("[copykat] getwd()=", getwd())
  message("[copykat] dir.exists('.')=", dir.exists("."))
  message("[copykat] list.files('.')= ", paste(list.files(".", all.files=TRUE), collapse=", "))
  quit(save="no", status = 2)
})

## -- labels（含 is_ref 审计列） --
anno_df <- NULL
if (!inherits(fit, "try-error") && !is.null(fit$prediction)) {
  anno_df <- fit$prediction; colnames(anno_df) <- c("barcode","label")
} else {
  anno_df <- data.frame(barcode = colnames(expr), label = "unknown", stringsAsFactors = FALSE)
}
anno_df$is_ref <- anno_df$barcode %in% ref_matched
write.csv(anno_df, "copykat_label.csv", row.names = FALSE)

file.create("COPYKAT.OK")
message("[copykat] done.")