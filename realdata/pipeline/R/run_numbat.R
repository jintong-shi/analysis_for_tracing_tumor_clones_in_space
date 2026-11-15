#!/usr/bin/env Rscript

sargs <- commandArgs(trailingOnly = TRUE)
getArg <- function(key, default=NULL) {
  hit <- grep(paste0("^", key, "="), sargs, value=TRUE)
  if (length(hit)==0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

visium     <- getArg("--visium_folder")
barcodes   <- getArg("--barcodes_file")
ref_barc   <- getArg("--ref_barcodes_file")
allele_tsv <- getArg("--allele_tsv_gz")
outdir     <- getArg("--outdir")
genome     <- getArg("--genome", "hg38")
ncores     <- as.integer(getArg("--ncores","4"))
max_ent    <- as.numeric(getArg("--max_entropy","0.8"))
min_LLR    <- as.numeric(getArg("--min_LLR","5"))
t_val      <- as.numeric(getArg("--t","1e-5"))        # ★ 新增
plot_flag  <- tolower(getArg("--plot","false")) %in% c("1","true","t","yes","y")
min_cells  <- as.integer(getArg("--min_cells","5")) 

if (any(is.null(c(visium,barcodes,ref_barc,allele_tsv,outdir)))) {
  stop("Missing required args. Need --visium_folder, --barcodes_file, --ref_barcodes_file, --allele_tsv_gz, --outdir")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(data.table)
  library(dplyr)
  library(numbat)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 读条形码
keep <- readLines(barcodes, warn = FALSE)
keep <- keep[nzchar(keep)]
ref  <- readLines(ref_barc, warn = FALSE)
ref  <- ref[nzchar(ref)]
if (length(keep) == 0) stop("keep_barcodes_file has no barcodes.")
if (length(ref)  == 0) stop("ref_barcodes_file has no barcodes (Numbat with-ref requires ref).")

# 读 Visium（10x）
seu <- Load10X_Spatial(visium)  # 默认使用 counts
cells <- colnames(seu)
keep <- intersect(keep, cells)
ref  <- intersect(ref,  cells)
if (length(keep) == 0) stop("None of keep barcodes are in Visium matrix.")
if (length(ref)  == 0) stop("None of ref barcodes are in Visium matrix.")

ref_seu   <- subset(seu, cells = ref)
keep_seu  <- subset(seu, cells = keep)
ref_counts  <- GetAssayData(ref_seu,  layer="counts")
keep_counts <- GetAssayData(keep_seu, layer="counts")

cell_annot <- data.frame(cell = colnames(ref_counts), group = "ref")
ref_internal <- aggregate_counts(ref_counts, cell_annot)

# 读等位数据
df_allele <- read_tsv(allele_tsv, show_col_types = FALSE)

message("[Numbat] running ...")
out <- tryCatch({
  run_numbat(
    keep_counts, ref_internal, df_allele,
    genome = genome,
    max_entropy = max_ent,
    t = t_val,              
    ncores = ncores,
    plot = plot_flag,
    min_cells = min_cells,
    min_LLR = min_LLR,
    out_dir = outdir
  )
}, error = function(e) {
  message("[Numbat] ERROR: ", e$message)
  quit(save="no", status=3)
})

# 判定是否“无 CNV”
bulk_file <- file.path(outdir, "bulk_clones_final.tsv.gz")
has_cnv <- FALSE
if (file.exists(bulk_file) && file.size(bulk_file) > 0) {
  bf <- tryCatch(fread(bulk_file), error=function(e) NULL)
  if (!is.null(bf) && "cnv_state" %in% names(bf)) {
    st <- bf$cnv_state[!is.na(bf$cnv_state)]
    has_cnv <- any(st %in% c("amp","del","bamp","bdel"))
  }
}

if (!has_cnv) {
  writeLines("No CNV remains after filtering by LLR in pseudobulks (or no non-neutral calls).",
             file.path(outdir, "NUMBAT.NO_CNV"))
  message("[Numbat] No CNV. Wrote NUMBAT.NO_CNV")
  quit(save="no", status=0)
}

# 正常导出：clone_df / label_df / bulk_for_wes.tsv
nb <- Numbat$new(out_dir = outdir)

clone_df <- nb$clone_post[, c("cell","clone_opt")]
fwrite(clone_df, file.path(outdir, "clone_df.tsv"), sep="\t")

label_df <- nb$clone_post[, c("cell","compartment_opt")]
fwrite(label_df, file.path(outdir, "label_df.tsv"), sep="\t")

dt <- as.data.table(nb$bulk_clones)
bulk_wes <- dt[, .(
  chr   = paste0("chr", as.character(CHROM)),
  start = seg_start,
  end   = seg_end,
  logfc = logFC,
  cnv_state_post,
  snp_id,
  gene
)]
fwrite(bulk_wes, file.path(outdir, "numbat_bulk_for_wes.tsv"), sep="\t", na="NA", quote=FALSE)

message("[Numbat] DONE and exported.")