#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(readr)
  library(magick)
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(FNN)
  library(igraph)
  library(TUSCAN)
  # 不加载 R 的 keras 包，使用 TUSCAN 内部的 reticulate/keras 后端
})

option_list <- list(
  make_option(c("--visium_folder"),     type = "character", help = "Path to Visium folder (root or outs/)"),
  make_option(c("--output_dir"),        type = "character", help = "Output dir for TUSCAN results"),
  make_option(c("--threads"),           type = "integer",   default = 8),
  make_option(c("--keep_barcodes_file"),type = "character", default = NULL),
  make_option(c("--seurat_rds"),        type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$visium_folder) || is.null(opt$output_dir)) {
  stop("--visium_folder and --output_dir are required")
}
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$output_dir)

get_counts_from_visium <- function(visium_folder) {
  mtx_dirs <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix"),
    file.path(visium_folder, "filtered_feature_bc_matrix"),
    visium_folder
  )
  has_mtx <- function(d) file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))
  mtx_dir <- mtx_dirs[file.exists(mtx_dirs) & vapply(mtx_dirs, has_mtx, logical(1))]
  if (length(mtx_dir) > 0) {
    message("[tuscan] Using 10x MTX directory: ", mtx_dir[1])
    return(Seurat::Read10X(mtx_dir[1]))
  }

  h5_files <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_files <- h5_files[file.exists(h5_files)]
  if (length(h5_files) > 0) {
    message("[tuscan] Using 10x H5 file: ", h5_files[1])
    return(Seurat::Read10X_h5(h5_files[1]))
  }

  stop("Neither 10x MTX directory nor filtered_feature_bc_matrix.h5 found under visium_folder: ", visium_folder)
}

read_tissue_positions <- function(visium_folder) {
  spatial_dirs <- c(
    file.path(visium_folder, "outs", "spatial"),
    file.path(visium_folder, "spatial"),
    file.path(dirname(visium_folder), "spatial")
  )
  spatial_dirs <- spatial_dirs[dir.exists(spatial_dirs)]
  if (length(spatial_dirs) == 0) stop("[tuscan] spatial/ folder not found under visium_folder or parent")

  sp <- spatial_dirs[1]
  list_path   <- file.path(sp, "tissue_positions_list.csv")
  header_path <- file.path(sp, "tissue_positions.csv")

  if (file.exists(header_path)) {
    df <- read.csv(header_path, header = TRUE, check.names = FALSE)
    # 规范列名
    col_map <- c(
      barcode = "barcode",
      in_tissue = "in_tissue",
      array_row = "array_row",
      array_col = "array_col",
      pxl_row_in_fullres = "pxl_row",
      pxl_col_in_fullres = "pxl_col",
      pxl_row = "pxl_row",
      pxl_col = "pxl_col"
    )
    names(df) <- ifelse(names(df) %in% names(col_map), col_map[names(df)], names(df))
    need <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
    missing <- setdiff(need, names(df))
    if (length(missing)) stop("[tuscan] tissue_positions.csv missing columns: ", paste(missing, collapse = ","))
    df <- df[, need]
    message("[tuscan] Using tissue positions file: ", header_path)
    return(df)
  }

  if (file.exists(list_path)) {
    df <- read.csv(list_path, header = FALSE)
    if (ncol(df) != 6) stop("[tuscan] tissue_positions_list.csv must have 6 columns, got: ", ncol(df))
    colnames(df) <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
    message("[tuscan] Using tissue positions file: ", list_path)
    return(df)
  }

  stop("[tuscan] Neither tissue_positions.csv nor tissue_positions_list.csv found under: ", sp)
}

# ---------- 读数据 ----------
# counts
counts <- get_counts_from_visium(opt$visium_folder)

# tissue positions
# tp_candidates <- c(
#   file.path(opt$visium_folder, "spatial", "tissue_positions.csv"),
#   file.path(opt$visium_folder, "spatial", "tissue_positions_list.csv"),
#   file.path(opt$visium_folder, "outs", "spatial", "tissue_positions.csv"),
#   file.path(opt$visium_folder, "outs", "spatial", "tissue_positions_list.csv")
# )
#tp_file <- tp_candidates[file.exists(tp_candidates)][1]
#if (is.na(tp_file)) stop("tissue_positions(.csv/_list.csv) not found")
#loc_raw <- suppressMessages(readr::read_csv(tp_file, col_names = TRUE, show_col_types = FALSE))
loc_raw <- read_tissue_positions(opt$visium_folder)

# read_tissue_positions() 已经把列名标准化为：barcode, pxl_row, pxl_col
# 这里直接用它们就行
required_cols <- c("barcode", "pxl_row", "pxl_col")
missing <- setdiff(required_cols, colnames(loc_raw))
if (length(missing)) {
  stop("[tuscan] tissue positions missing columns after standardization: ", paste(missing, collapse = ", "))
}

loc <- data.frame(
  x = loc_raw$pxl_col,
  y = loc_raw$pxl_row,
  row.names = loc_raw$barcode
)

# image & scale  —— 先找 hires，若无则用 lowres
img_candidates_hires <- c(
  file.path(opt$visium_folder, "spatial", "tissue_hires_image.png"),
  file.path(opt$visium_folder, "outs", "spatial", "tissue_hires_image.png")
)
img_candidates_lowres <- c(
  file.path(opt$visium_folder, "spatial", "tissue_lowres_image.png"),
  file.path(opt$visium_folder, "outs", "spatial", "tissue_lowres_image.png")
)

img_hires <- img_candidates_hires[file.exists(img_candidates_hires)][1]
img_low   <- img_candidates_lowres[file.exists(img_candidates_lowres)][1]

sf_candidates <- c(
  file.path(opt$visium_folder, "spatial", "scalefactors_json.json"),
  file.path(opt$visium_folder, "outs", "spatial", "scalefactors_json.json")
)
sf_file <- sf_candidates[file.exists(sf_candidates)][1]
if (is.na(sf_file)) stop("scalefactors_json.json not found")
scales <- jsonlite::fromJSON(sf_file)

if (!is.na(img_hires)) {
  img <- magick::image_read(img_hires)
  scale_val <- if (!is.null(scales$tissue_hires_scalef)) scales$tissue_hires_scalef else 1.0
  message("[tuscan] Using hires image: ", img_hires, "  (scale=", scale_val, ")")
} else if (!is.na(img_low)) {
  img <- magick::image_read(img_low)
  scale_val <- if (!is.null(scales$tissue_lowres_scalef)) scales$tissue_lowres_scalef else 1.0
  message("[tuscan] Using lowres image: ", img_low, "  (scale=", scale_val, ")")
} else {
  stop("Neither tissue_hires_image.png nor tissue_lowres_image.png found under spatial/")
}

# loc_* 是 fullres 坐标，需乘当前图像的缩放因子
loc$x <- loc$x * scale_val
loc$y <- loc$y * scale_val

# 对齐条形码
common <- intersect(colnames(counts), rownames(loc))
counts <- counts[, common, drop = FALSE]
loc    <- loc[common, , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(loc)))

# ---------- 只用 keep_barcodes 做子集（唯一的筛选来源） ----------
if (!is.null(opt$keep_barcodes_file) && file.exists(opt$keep_barcodes_file)) {
  keep <- readr::read_lines(opt$keep_barcodes_file)
  keep <- intersect(keep, colnames(counts))
  if (length(keep) == 0) stop("After applying keep_barcodes_file, no barcodes remain.")
  counts <- counts[, keep, drop = FALSE]
  loc    <- loc[keep, , drop = FALSE]
  stopifnot(all(colnames(counts) == rownames(loc)))
}

# ---------- TUSCAN 主流程 ----------
set.seed(1234)
obj <- createTUSCANObject(counts = counts, location = loc, img = img, project = "TUSCAN", r = 49)
obj <- reduceDim(object = obj, nhvg = 2000, d = 20, epochs = 50, batch_size = 64, seed = 1234)

# KNN 图 + 最优分辨率
compute_graph <- function(features, k = 100) {
  knn.norm <- FNN::get.knn(as.matrix(features), k = k)
  knn.norm <- data.frame(from = rep(1:nrow(knn.norm$nn.index), k),
                         to = as.vector(knn.norm$nn.index),
                         weight = 1 / (1 + as.vector(knn.norm$nn.dist)))
  nw.norm <- igraph::graph_from_data_frame(knn.norm, directed = FALSE)
  igraph::simplify(nw.norm)
}
find_best_resolution_modularity <- function(graph, res_values = seq(0.1, 1.5, by = 0.1)) {
  best_res <- res_values[1]; best_mod <- -Inf
  for (res in res_values) {
    cl <- igraph::cluster_louvain(graph, resolution = res)
    mod <- igraph::modularity(graph, cl$membership)
    if (mod > best_mod) { best_mod <- mod; best_res <- res }
  }
  best_res
}

features <- cbind(obj@autoencoder, obj@image[, c("RS", "GS", "BS")])
graph <- compute_graph(features)
opt_res <- find_best_resolution_modularity(graph)
obj <- PreCluster(object = obj, res = opt_res, method = "louvain", seed = 1234)

ref_group <- FindNormalCluster(obj)
p <- plotSpatial(obj = obj, class = "cluster", size = 1) + ggplot2::ggtitle(paste0("ref_group=", ref_group))
ggsave(filename = file.path(opt$output_dir, "PreCluster.pdf"), plot = p, width = 7, height = 6.5)

cnv_obj <- runCNV(object = obj,
                  ref_group_names = ref_group,
                  min_ave_counts_per_gene = 0.01,
                  min_spots_per_gene = 3,
                  chr_exclude = c("chrX", "chrY", "chrM"),
                  normalization = FALSE)

plotCNVheatmap(object = cnv_obj, title = "heatmap", output_format = "png")

tumor <- TumorCluster(object = cnv_obj, seed = 1234)
p1 <- plotSpatial(obj = tumor, class = "tumor", size = 1)
ggsave(filename = file.path(opt$output_dir, "TumorCluster.pdf"), plot = p1, width = 7, height = 6.5)

# 保存对象
saveRDS(obj,     file = file.path(opt$output_dir, "obj.rds.gz"),     compress = "gzip")
saveRDS(cnv_obj, file = file.path(opt$output_dir, "cnv_obj.rds.gz"), compress = "gzip")
saveRDS(tumor,   file = file.path(opt$output_dir, "tumor_obj.rds.gz"),compress = "gzip")

# 导出标签（便于外部整合）
try({
  labs  <- tumor@tumor_cluster
  labs2 <- ifelse(labs %in% c("Tumor","tumor","aneuploid"), "aneuploid",
                  ifelse(labs %in% c("Normal","normal","diploid"), "diploid", NA))
  anno_df <- data.frame(barcode = names(labs2), label = labs2)
  write.csv(anno_df, file = file.path(opt$output_dir, "TUSCAN_label.csv"), row.names = FALSE)
}, silent = TRUE)

# 可选：画个 Seurat 的 SpatialDimPlot
seu <- NULL
if (!is.null(opt$seurat_rds) && file.exists(opt$seurat_rds)) {
  seu <- try(readRDS(opt$seurat_rds), silent = TRUE)
} else if (!is.null(opt$visium_folder)) {
  data.dir <- opt$visium_folder
  if (file.exists(file.path(opt$visium_folder, "outs", "spatial"))) {
    data.dir <- file.path(opt$visium_folder, "outs")
  }
  seu <- try(Seurat::Load10X_Spatial(data.dir = data.dir), silent = TRUE)
}
if (!is.null(seu) && !inherits(seu, "try-error")) {
  try({
    lab <- setNames(anno_df$label, anno_df$barcode)
    seu$aneu_label <- lab[colnames(seu)]
    p <- try(Seurat::SpatialDimPlot(seu, group.by = "aneu_label",
                                    cols = c(aneuploid = "#d62728", diploid = "#2ca02c")),
             silent = TRUE)
    if (!inherits(p, "try-error")) {
      ggplot2::ggsave(filename = file.path(opt$output_dir, "SpatialDimplot_TUSCAN.pdf"),
                      plot = p, width = 7, height = 6.5)
    }
  }, silent = TRUE)
}

message("TUSCAN done. Output: ", opt$output_dir)
# 正确的哨兵文件名
file.create(file.path(opt$output_dir, "TUSCAN.OK"))
flush.console()