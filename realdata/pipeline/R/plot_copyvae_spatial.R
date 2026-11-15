#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--visium_folder"),      type="character", help="Path to Visium folder (root or outs)"),
  make_option(c("--output_dir"),         type="character", help="copyVAE output dir (contains pred_label.txt, barcodes_order.txt)"),
  make_option(c("--keep_barcodes_file"), type="character", help="Text file of kept barcodes (one per line)"),
  make_option(c("--plot_file"),          type="character", default="SpatialDimplot_copyvae.pdf",
              help="Output PDF filename [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$visium_folder), !is.null(opt$output_dir), !is.null(opt$keep_barcodes_file))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$output_dir)

# --- 定位 10x HDF5 矩阵 ---
data_dir <- opt$visium_folder
if (!file.exists(file.path(data_dir, "filtered_feature_bc_matrix.h5"))) {
  alt <- file.path(opt$visium_folder, "outs")
  if (file.exists(file.path(alt, "filtered_feature_bc_matrix.h5"))) data_dir <- alt
}
if (!file.exists(file.path(data_dir, "filtered_feature_bc_matrix.h5")))
  stop("filtered_feature_bc_matrix.h5 not found under visium folder or visium/outs")

# --- 读 Visium + 图片 ---
seu_full <- Seurat::Load10X_Spatial(
  data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5", filter.matrix = TRUE
)
barcodes_all <- colnames(seu_full)

# --- 读 copyVAE 结果（对齐 barcodes_order.txt） ---
order_file <- file.path(opt$output_dir, "barcodes_order.txt")
if (!file.exists(order_file)) stop("barcodes_order.txt not found in output_dir; rerun copyVAE prepare")
barcodes_order <- readr::read_lines(order_file)

pred_path <- file.path(opt$output_dir, "pred_label.txt")
if (!file.exists(pred_path)) stop("pred_label.txt not found in output_dir")
pred <- read.table(pred_path, header = FALSE)
colnames(pred) <- "cluster"

# 长度对齐到 barcodes_order
if (nrow(pred) != length(barcodes_order)) {
  warning(sprintf("pred rows (%d) != barcodes_order (%d); truncation applied",
                  nrow(pred), length(barcodes_order)))
  n <- min(nrow(pred), length(barcodes_order))
  pred <- pred[seq_len(n), , drop = FALSE]
  barcodes_order <- barcodes_order[seq_len(n)]
}
rownames(pred) <- barcodes_order

# 正常类索引（可选）
ncell_file <- file.path(opt$output_dir, "ncell_index.csv")
normal_idx <- 0
if (file.exists(ncell_file)) {
  normal_idx <- tryCatch({ as.integer(read.csv(ncell_file)[1,1]) }, error = function(e) 0)
}

# --- 只用 keep_barcodes_file 来筛选 ---
keep_all <- unique(readr::read_lines(opt$keep_barcodes_file))
keep_all <- keep_all[nchar(keep_all) > 0]

# 只保留：既在 barcodes_order（copyVAE 参与的细胞），也在原 Visium 中存在的
keep <- intersect(keep_all, intersect(barcodes_order, barcodes_all))
if (length(keep) == 0) stop("No overlapping barcodes among keep list, copyVAE order, and Visium data")

# 从 pred 里切换到 keep 顺序
pred_sub <- pred[keep, , drop = FALSE]
pred_sub$copyvae <- ifelse(pred_sub$cluster == normal_idx, "diploid", "aneuploid")

# --- 子集并加上元数据 ---
seu <- subset(seu_full, cells = keep)
seu <- AddMetaData(seu, metadata = pred_sub["copyvae"])  # 传 data.frame

p <- SpatialDimPlot(seu, group.by = "copyvae", pt.size.factor = 2) +
  ggplot2::ggtitle("copyVAE: diploid vs aneuploid")
ggsave(filename = file.path(opt$output_dir, opt$plot_file), plot = p, width = 7, height = 6.5)

invisible(TRUE)