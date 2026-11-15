#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(readr)
})

opt_list <- list(
  make_option(c("--visium_folder"), type="character", help="Visium 根目录或 outs/ 目录"),
  make_option(c("--out_keep"),      type="character", help="输出 keep_barcodes.txt")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$visium_folder) || is.null(opt$out_keep)) {
  stop("--visium_folder 和 --out_keep 均为必填")
}

# 找 10x 数据（优先 MTX 目录；否则 H5）
candidate_dirs <- c(
  file.path(opt$visium_folder, "outs", "filtered_feature_bc_matrix"),
  file.path(opt$visium_folder, "filtered_feature_bc_matrix"),
  opt$visium_folder
)
has_mtx <- function(d) file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))
mtx_dirs <- candidate_dirs[dir.exists(candidate_dirs) & vapply(candidate_dirs, has_mtx, logical(1))]

counts <- NULL
if (length(mtx_dirs) > 0) {
  message("[allspots] Using 10x MTX directory: ", mtx_dirs[1])
  counts <- Seurat::Read10X(mtx_dirs[1])
} else {
  h5_candidates <- c(
    file.path(opt$visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(opt$visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(opt$visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_candidates <- h5_candidates[file.exists(h5_candidates)]
  if (length(h5_candidates) == 0)
    stop("Neither MTX dir nor filtered_feature_bc_matrix.h5 found under visium_folder: ", opt$visium_folder)
  message("[allspots] Using 10x H5 file: ", h5_candidates[1])
  counts <- Seurat::Read10X_h5(h5_candidates[1])
}

barcodes <- colnames(counts)
if (length(barcodes) == 0) stop("No barcodes found in Visium data.")

dir.create(dirname(opt$out_keep), showWarnings = FALSE, recursive = TRUE)
readr::write_lines(barcodes, opt$out_keep)
message("[allspots] Wrote ", length(barcodes), " barcodes -> ", opt$out_keep)
quit(save = "no", status = 0)