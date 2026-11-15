#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SCEVAN)
  library(readr)
})

# SCEVAN's plotTSNE() can abort with "data.frame: arguments imply differing
# number of rows" when the subclone tSNE has fewer points than predicted labels.
# It is called one line *before* the subclone column is attached to classDf, so
# the failure otherwise discards already-computed subclone assignments. Guard it
# so pipelineCNA still returns the labelled result.
if ("plotTSNE" %in% ls(getNamespace("SCEVAN"))) {
  .orig_plotTSNE <- get("plotTSNE", envir = getNamespace("SCEVAN"))
  assignInNamespace("plotTSNE", function(...) tryCatch(
    .orig_plotTSNE(...),
    error = function(e) {
      message("[scevan] plotTSNE skipped: ", conditionMessage(e)); invisible(NULL)
    }
  ), ns = "SCEVAN")
}

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
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

## -- 统一创建并切到输出目录（绝对路径） --
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
outdir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
message("[scevan] outdir=", outdir)
if (!dir.exists(outdir)) stop("Output dir not found after dir.create: ", outdir)
old_wd <- getwd()
setwd(outdir)
on.exit(setwd(old_wd), add = TRUE)

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
    message("[SCEVAN] Using 10x MTX directory: ", mtx_dir[1])
    return(Seurat::Read10X(mtx_dir[1]))
  }
  h5_files <- c(
    file.path(visium_folder, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.h5"),
    file.path(visium_folder, "filtered_feature_bc_matrix.hdf5")
  )
  h5_files <- h5_files[file.exists(h5_files)]
  if (length(h5_files) > 0) {
    message("[SCEVAN] Using 10x H5 file: ", h5_files[1])
    return(Seurat::Read10X_h5(h5_files[1]))
  }
  stop("Neither 10x MTX dir nor filtered_feature_bc_matrix.h5 found: ", visium_folder)
}

trim_nonempty <- function(x) unique(x[nchar(x) > 0])

# --- load counts ---
counts <- get_counts_from_visium(opt$visium_folder)

# --- read keep/ref (raw) ---
keep_raw <- trim_nonempty(readr::read_lines(opt$keep_barcodes_file))
ref_raw  <- character(0)
if (!is.null(opt$ref_barcodes_file) && nzchar(opt$ref_barcodes_file) && file.exists(opt$ref_barcodes_file)) {
  ref_raw <- trim_nonempty(readr::read_lines(opt$ref_barcodes_file))
}

# --- match to matrix columns ---
keep_matched <- intersect(keep_raw, colnames(counts))
ref_matched  <- intersect(ref_raw,  colnames(counts))

# --- AUDIT ---
audit_dir <- opt$output_dir
writeLines(c(
  sprintf("keep_path: %s", normalizePath(opt$keep_barcodes_file, winslash="/", mustWork=FALSE)),
  sprintf(" ref_path: %s", if (is.null(opt$ref_barcodes_file)) "(none)" else normalizePath(opt$ref_barcodes_file, winslash="/", mustWork=FALSE)),
  sprintf("keep_raw: %d", length(keep_raw)),
  sprintf(" ref_raw: %d", length(ref_raw)),
  sprintf("keep_matched: %d", length(keep_matched)),
  sprintf(" ref_matched: %d", length(ref_matched))
), file.path(audit_dir, "REF_AUDIT.txt"))
writeLines(keep_matched, file.path(audit_dir, "keep_used.txt"))
writeLines(ref_matched,  file.path(audit_dir, "ref_used.txt"))
message(sprintf("[SCEVAN][audit] keep_raw=%d keep_matched=%d; ref_raw=%d ref_matched=%d",
                length(keep_raw), length(keep_matched), length(ref_raw), length(ref_matched)))

# --- select columns = union(ref, keep) preserving original order ---
sel <- union(ref_matched, keep_matched)
sel <- colnames(counts)[colnames(counts) %in% sel]
if (length(sel) == 0) stop("No cells selected after union(ref, keep).")
expr <- counts[, sel, drop = FALSE]

# --- run SCEVAN (with/without ref) ---
set.seed(1234)
res <- try(
  if (length(ref_matched) > 0) {
    pipelineCNA(expr, par_cores = max(1, opt$threads - 1),
                plotTree = TRUE, organism = 'human', norm_cell = ref_matched)
  } else {
    pipelineCNA(expr, par_cores = max(1, opt$threads - 1),
                plotTree = TRUE, organism = 'human')
  },
  silent = TRUE
)
saveRDS(res, file = file.path(opt$output_dir, "SCEVAN_results.rds"))

# --- label CSV (barcode, label, is_ref) ---
anno_df <- NULL
if (!inherits(res, "try-error")) {
  df <- try(as.data.frame(res), silent=TRUE)
  if (!inherits(df, "try-error") && nrow(df) > 0) {
    if (!"class" %in% colnames(df)) colnames(df)[1] <- "class"
    df$label <- ifelse(df$class %in% c("tumor","Tumor","aneuploid"), "aneuploid",
                       ifelse(df$class %in% c("normal","Normal","diploid"), "diploid", "unknown"))
    anno_df <- data.frame(barcode = rownames(df), label = df$label, stringsAsFactors=FALSE)
  }
}
if (is.null(anno_df)) {
  anno_df <- data.frame(barcode = colnames(expr), label = "unknown", stringsAsFactors=FALSE)
}
anno_df$is_ref <- anno_df$barcode %in% ref_matched
write.csv(anno_df, file.path(opt$output_dir, "SCEVAN_label.csv"), row.names = FALSE)

file.create(file.path(opt$output_dir, "SCEVAN.OK"))
message("[SCEVAN] done.")