#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(readr)
  library(infercnv)
})

opt_list <- list(
  make_option(c("--visium_folder"), type="character"),
  make_option(c("--output_dir"), type="character"),
  make_option(c("--keep_barcodes_file"), type="character"),
  make_option(c("--ref_barcodes_file"), type="character", default=NULL),
  make_option(c("--gene_pos_file"), type="character"),
  make_option(c("--threads"), type="integer", default=8),
  make_option(c("--hmm"), type="character", default=TRUE) 
)
opt <- parse_args(OptionParser(option_list=opt_list))

for (k in c("visium_folder","output_dir","keep_barcodes_file","gene_pos_file"))
  if (is.null(opt[[k]])) stop("--", k, " is required")
if (!file.exists(opt$gene_pos_file)) stop("gene_pos_file not found: ", opt$gene_pos_file)
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 决定是否跑 HMM（默认 FALSE）----
str2bool <- function(x) tolower(x) %in% c("true","t","1","yes","y")
if (is.null(opt$hmm)) {
  do_hmm <- FALSE
} else {
  do_hmm <- str2bool(opt$hmm)
}

# ---------- IO helpers ----------
trim_nonempty <- function(x) unique(x[nchar(x) > 0])

get_counts_from_visium <- function(visium_folder) {
  mtx_dirs <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix"),
    file.path(visium_folder, "filtered_feature_bc_matrix"),
    visium_folder
  )
  has_mtx <- function(d) {
    file.exists(file.path(d, "matrix.mtx")) || file.exists(file.path(d, "matrix.mtx.gz"))
  }
  mtx_dir <- mtx_dirs[file.exists(mtx_dirs) & vapply(mtx_dirs, has_mtx, logical(1))]
  if (length(mtx_dir) > 0) {
    message("[infercnv] Using 10x MTX directory: ", mtx_dir[1])
    return(Seurat::Read10X(mtx_dir[1]))
  }

  h5_files <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_files <- h5_files[file.exists(h5_files)]
  if (length(h5_files) > 0) {
    message("[infercnv] Using 10x H5 file: ", h5_files[1])
    return(Seurat::Read10X_h5(h5_files[1]))
  }

  stop("Neither 10x MTX directory nor filtered_feature_bc_matrix.h5 found under: ", visium_folder)
}

# ---------- barcode harmonizer (handles '-1') ----------
harmonize_ids <- function(colv, ids) {
  # colv: colnames(counts), ids: input keep/ref vectors (may be empty)
  ids <- trim_nonempty(ids)
  if (!length(ids)) return(character(0))

  # 1) direct match
  hit <- intersect(ids, colv)
  if (length(hit)) return(hit)

  # 2) try add "-1"
  hit <- intersect(paste0(ids, "-1"), colv)
  if (length(hit)) return(hit)

  # 3) try strip "-1" on both sides, then map back to original colv
  strip1 <- sub("-1$", "", ids)
  strip2 <- sub("-1$", "", colv)
  hit_stripped <- intersect(strip1, strip2)
  if (length(hit_stripped)) {
    mask <- strip2 %in% hit_stripped
    return(colv[mask])
  }

  character(0)
}

# ---------- read counts ----------
counts <- get_counts_from_visium(opt$visium_folder)
all_cols <- colnames(counts)

# ---------- read keep/ref ----------
keep_raw <- trim_nonempty(readr::read_lines(opt$keep_barcodes_file))
ref_raw  <- character(0)
if (!is.null(opt$ref_barcodes_file) && nzchar(opt$ref_barcodes_file) && file.exists(opt$ref_barcodes_file)) {
  ref_raw <- trim_nonempty(readr::read_lines(opt$ref_barcodes_file))
}

# ---------- harmonize against counts colnames ----------
keep_h <- harmonize_ids(all_cols, keep_raw)
ref_h  <- harmonize_ids(all_cols, ref_raw)

# 日志打印
msg <- function(...) cat(sprintf(...), "\n")
msg("[infercnv] keep: raw=%d, matched=%d", length(keep_raw), length(keep_h))
msg("[infercnv]  ref: raw=%d, matched=%d", length(ref_raw),  length(ref_h))
msg("[infercnv] hmm:",ifelse(do_hmm,"TRUE","FALSE"))


# ---------- subset by union(ref, keep), preserving original order ----------
sel <- union(ref_h, keep_h)
sel <- all_cols[all_cols %in% sel]  # keep original order
if (length(sel) == 0) {
  stop("No cells selected after harmonization. Check keep/ref barcodes vs matrix columns.")
}
counts <- counts[, sel, drop = FALSE]
cells  <- colnames(counts)

# 以最终 cells 再裁剪 ref_h，确保 ref ⊆ cells
ref_final <- intersect(ref_h, cells)

# ---------- annotations ----------
annot_file <- file.path(opt$output_dir, "annotations.tsv")
if (length(ref_final) > 0) {
  is_ref <- cells %in% ref_final
  annot_df <- data.frame(cell = cells,
                         group = ifelse(is_ref, "Normal", "Tumor"),
                         stringsAsFactors = FALSE)
  ref_names <- c("Normal")
} else {
  annot_df <- data.frame(cell = cells,
                         group = rep("Tumor", length(cells)),
                         stringsAsFactors = FALSE)
  ref_names <- NULL
}
write.table(annot_df, annot_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# 简要统计保存，方便你排查
write.table(as.data.frame(table(annot_df$group)),
            file.path(opt$output_dir, "_group_counts.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# ---------- run infercnv ----------
# Keep no-ref compatibility, but follow the requested parameter set when ref exists.
if (do_hmm && length(ref_final) == 0) {
  warning("HMM requested but no Normal reference barcodes matched; disabling HMM.")
  do_hmm <- FALSE
}

if (length(ref_final) == 0) {
  warning("No Normal reference barcodes matched; running infercnv without ref_group_names.")
}

obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = annot_file,
  gene_order_file   = opt$gene_pos_file,
  delim             = "\t",
  ref_group_names   = ref_names
)

message(paste0("Running inferCNV with HMM = ", do_hmm))
res <- infercnv::run(
  obj,
  cutoff = 0.1,
  out_dir = opt$output_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = do_hmm,
  HMM_type = "i3",
  num_threads = opt$threads,
  HMM_report_by = "cell",
  no_prelim_plot = TRUE,
  no_plot = TRUE
)

file.create(file.path(opt$output_dir, "infercnv.OK"))
quit(save = "no", status = 0)
