#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(readr)
  library(dplyr)
  library(CONICSmat)
  library(data.table)
  library(purrr)
  library(Matrix)
  library(matrixStats)
  library(zoo)
  library(pdist)
})

opt_list <- list(
  make_option(c("--visium_folder"),      type="character"),
  make_option(c("--output_dir"),         type="character"),
  make_option(c("--keep_barcodes_file"), type="character"),
  make_option(c("--ref_barcodes_file"),  type="character"),
  make_option(c("--chrom_arm_file"),     type="character"),
  make_option(c("--gene_pos_file"),      type="character", default=""),
  make_option(c("--threads"),            type="integer",   default=8),
  make_option(c("--exp_thresh"),         type="double",   default=0.1),
  make_option(c("--windowsize"),         type="integer",   default=101),
  make_option(c("--threshold"),         type="double",   default=4)

)
ori_path = getwd()
opt <- parse_args(OptionParser(option_list=opt_list))
for (k in c("visium_folder","output_dir","keep_barcodes_file","chrom_arm_file")) {
  if (is.null(opt[[k]])) stop("--",k," is required")
}
if (!file.exists(opt$chrom_arm_file)) stop("chrom_arm_file not found: ", opt$chrom_arm_file)
if (nzchar(opt$gene_pos_file) && !file.exists(opt$gene_pos_file)) {
  stop("gene_pos_file not found: ", opt$gene_pos_file)
}
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$output_dir)

get_counts_from_visium <- function(visium_folder) {
  # 候选的 10x mtx 目录
  mtx_dirs <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix"),
    file.path(visium_folder, "filtered_feature_bc_matrix"),
    visium_folder  # 允许直接传入矩阵目录本身
  )
  has_mtx <- function(d) {
    file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))
  }
  mtx_dir <- mtx_dirs[file.exists(mtx_dirs) & vapply(mtx_dirs, has_mtx, logical(1))]
  if (length(mtx_dir) > 0) {
    message("[scevan] Using 10x MTX directory: ", mtx_dir[1])
    return(Seurat::Read10X(mtx_dir[1]))
  }

  # 候选的 h5 文件
  h5_files <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_files <- h5_files[file.exists(h5_files)]
  if (length(h5_files) > 0) {
    message("[scevan] Using 10x H5 file: ", h5_files[1])
    return(Seurat::Read10X_h5(h5_files[1]))
  }

  stop("Neither 10x MTX directory nor filtered_feature_bc_matrix.h5 found under visium_folder: ", visium_folder)
}

# --- Read Visium matrix (genes x spots) ---
counts <- get_counts_from_visium(opt$visium_folder)


# --- keep subsetting (empty file => use all) ---
keep <- character(0)
if (file.exists(opt$keep_barcodes_file)) {
  keep <- readr::read_lines(opt$keep_barcodes_file)
  keep <- intersect(keep, colnames(counts))
}
ref <- character(0)
if (file.exists(opt$ref_barcodes_file)) {
  ref <- readr::read_lines(opt$ref_barcodes_file)
  ref <- intersect(ref, colnames(counts))
}

if (length(keep) > 0) {
  counts <- counts[, unique(c(keep,ref)), drop = FALSE]
}
if (ncol(counts) < 50) stop("Too few spots after keep filtering (", ncol(counts), ").")

# --- regions (chromosome arms) ---
regions <- read.table(opt$chrom_arm_file, sep = "\t", row.names = 1, header = TRUE)

# --- CPM & log2 ---
suva_expr = CONICSmat::normMat(as.matrix(counts))

load_gene_pos <- function(fp) {
  gp <- data.table::fread(fp, data.table = FALSE)
  if (ncol(gp) < 4) stop("Gene position file must have at least 4 columns.")

  nm <- tolower(names(gp))
  pick <- function(cands) {
    idx <- which(nm %in% cands)
    if (length(idx) > 0) return(idx[1])
    NA_integer_
  }

  gene_col  <- pick(c("gene", "hgnc_symbol", "symbol"))
  chr_col   <- pick(c("chr", "chromosome", "chromosome_name"))
  start_col <- pick(c("start", "start_position", "txstart"))
  end_col   <- pick(c("end", "end_position", "txend"))

  if (any(is.na(c(gene_col, chr_col, start_col, end_col)))) {
    gene_col <- 1
    chr_col <- 2
    start_col <- 3
    end_col <- 4
  }

  gene  <- gsub('"', "", as.character(gp[[gene_col]]))
  chr   <- gsub('"', "", as.character(gp[[chr_col]]))
  start <- suppressWarnings(as.numeric(gp[[start_col]]))
  end   <- suppressWarnings(as.numeric(gp[[end_col]]))

  chr <- sub("^chr", "", chr, ignore.case = TRUE)

  data.frame(
    ensembl_gene_id = gene,
    hgnc_symbol = gene,
    chromosome_name = chr,
    start_position = start,
    end_position = end,
    stringsAsFactors = FALSE
  )
}

load_gene_pos_for_conicsmat <- function(genes, local_fp = "") {
  if (nzchar(local_fp)) {
    message("-> Loading gene positions from local file...")
    return(load_gene_pos(local_fp))
  }

  message("-> Fetching gene positions (Ensembl)...")
  CONICSmat::getGenePositions(genes, ensembl_version = "https://www.ensembl.org")
}

gene_pos <- load_gene_pos_for_conicsmat(rownames(suva_expr), opt$gene_pos_file)

# --- filter & normalize ---
suva_expr <- CONICSmat::filterMatrix(suva_expr, gene_pos[, "hgnc_symbol"], minCells = 5)
normFactor <- CONICSmat::calcNormFactors(suva_expr)

# Drop arms whose gene count after filterMatrix is too small for plotChrEnichment.
# plotChrEnichment -> scale -> colMeans crashes when expmat has <2 genes on an arm
# (drops to a vector). Acrocentric p-arms 13p/14p/15p have 0 protein-coding genes.
.gp_kept <- gene_pos[gene_pos$hgnc_symbol %in% rownames(suva_expr), , drop = FALSE]
.gp_kept$chromosome_name <- sub("^chr", "", .gp_kept$chromosome_name, ignore.case = TRUE)
.genes_per_arm <- vapply(rownames(regions), function(.arm) {
  .chr <- as.character(regions[.arm, "Chrom"])
  .s   <- as.numeric(regions[.arm, "Start"])
  .e   <- as.numeric(regions[.arm, "End"])
  sum(.gp_kept$chromosome_name == .chr &
      .gp_kept$start_position >= .s &
      .gp_kept$start_position <= .e, na.rm = TRUE)
}, integer(1))
.dropped <- rownames(regions)[.genes_per_arm < 2]
if (length(.dropped) > 0) {
  message(sprintf("[CONICSmat] dropping %d empty arm(s) (<2 genes after minCells filter): %s",
                  length(.dropped), paste(.dropped, collapse = ",")))
  regions <- regions[.genes_per_arm >= 2, , drop = FALSE]
}

# --- run & plots ---
l <- CONICSmat::plotAll(suva_expr, normFactor, regions, gene_pos, "SUVA_CNVs")
saveRDS(suva_expr, file.path(opt$output_dir, "SUVA_expr.rds"))

pdf(file.path(opt$output_dir, "plotHistogram_k2.pdf"))
invisible(CONICSmat::plotHistogram(l, suva_expr, clusters = 2, zscoreThreshold = opt$threshold))
dev.off()

pdf(file.path(opt$output_dir, "plotChromosomeHeatmap.pdf"))
invisible(CONICSmat::plotChromosomeHeatmap(
  suva_expr,
  normal     = which(colnames(suva_expr)%in%ref),
  plotcells  = colnames(suva_expr),
  gene_pos   = gene_pos,
  windowsize = opt$windowsize,
  chr        = TRUE,
  expThresh  = opt$exp_thresh,
  thresh     = opt$threshold
))
dev.off()

file.create(file.path(opt$output_dir, "CONICSmat.OK"))
cat("CONICSmat run complete. Output in:", opt$output_dir, "\n")
invisible(TRUE)

# --- 参考：有 normal 用 normal 均值----
expThresh <- opt$exp_thresh
plotcells <- colnames(suva_expr)

ref_vec = rowMeans(suva_expr[, ref])

tumor = setdiff(plotcells, ref)
gexp = suva_expr[, c(ref, tumor)]

ref_vec = ref_vec[which(ref_vec > expThresh)]
gexp = gexp[names(ref_vec), ]
gexp = gexp[which(rowMeans(gexp) > expThresh), ]
ref_vec = ref_vec[rownames(gexp)]
gene_pos <- gene_pos %>%
  mutate(
    chromosome_name = gsub("^chr", "", chromosome_name, ignore.case = TRUE),
    chr_i   = suppressWarnings(as.integer(as.character(chromosome_name))),
    start_i = map_dbl(start_position, ~ suppressWarnings(as.double(.x[1])))
  ) %>%
  filter(chr_i %in% 1:24, is.finite(start_i)) %>%
  arrange(chr_i, start_i) %>%
  select(-chr_i, -start_i)
colnames(gene_pos)[2]=c('gene')
gexp = gexp[intersect(gene_pos$gene, rownames(gexp)), ]
ref_vec = ref_vec[rownames(gexp)]
rat = gexp - ref_vec
print("Generating distance Matrix")
d = apply(rat, 2, function(x) zoo::rollapply(x, mean, width = opt$windowsize, 
                                             align = "center"))
d = apply(d, 2, function(x) x - mean(x))
print(paste("Distance matrix ready with dimensions", nrow(d), 
            ncol(d)))

print("Clustering cells")
ref = 1:length(ref)
tumor = (length(ref) + 1):ncol(d)
dmat = as.dist(pdist(t(d[, tumor])))
hc = hclust(dmat)
cellOrder = c(ref, tumor[hc$order])
d = d[, cellOrder]

# 保存
gp <- dplyr::filter(gene_pos,gene%in%rownames(d))
axis_csv <- file.path(opt$output_dir, "CONICS_axis.csv")
saveRDS(d, file.path(opt$output_dir, "CONICS_track_matrix.rds"))
fwrite(gp, axis_csv)

# 也写一个 OK 标记
writeLines("OK", file.path(opt$output_dir, "CONICS_track.OK"))
setwd(ori_path)
