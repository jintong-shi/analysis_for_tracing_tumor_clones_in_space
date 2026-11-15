#!/usr/bin/env Rscript
# atera_20_cnvtools_prep.R
#
# Export the atera (Xenium) benign-reference cell set as a 10x MTX directory so
# the shared pipeline scripts snake_pipeline/R/run_conicsmat.R and run_copykat.R
# can be run on it unmodified.
#
# Counts come from the frozen inferCNV run's @count.data (raw integer counts,
# 7,951 genes x 87,989 cells) rather than the Xenium Seurat object, so CONICSmat
# and CopyKAT see exactly the cells and genes inferCNV saw:
#   ref  = @reference_grouped_cell_indices$normal   (2,799 cells)
#   keep = @observation_grouped_cell_indices$cancer (84,247 cells)
# Cells in neither group are dropped, as they were for inferCNV.
#
# Note: outputs/xenium_atera.rds is NOT usable here -- it holds CNV residuals for
# the cancer cells only (n_benign = 0), not raw counts.
#
# Env: stcnv-r-shared.
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/atera/code/atera_20_cnvtools_prep.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
})

## ---- config -----------------------------------------------------
ATERA    <- "<DATA_ROOT>/high_res_ST/atera"
obj_path <- file.path(ATERA, "processed_data/run4_benign_ref_strict_qc/run.final.infercnv_obj")
input_dir <- file.path(ATERA, "processed_data/cnvtools/atera/input")

## ---- 1. raw counts + group membership ---------------------------
message("reading inferCNV object (1.5 GB) ...")
o <- readRDS(obj_path)

counts <- o@count.data
ref    <- colnames(counts)[o@reference_grouped_cell_indices$normal]
keep   <- colnames(counts)[o@observation_grouped_cell_indices$cancer]

if (length(ref) == 0L || length(keep) == 0L)
  stop("expected both a 'normal' reference group and a 'cancer' observation group")

message(sprintf("ref (normal) = %d cells, keep (cancer) = %d cells, %d genes",
                length(ref), length(keep), nrow(counts)))

counts <- counts[, c(keep, ref), drop = FALSE]
counts <- as(as(Matrix(counts, sparse = TRUE), "CsparseMatrix"), "generalMatrix")
rm(o); invisible(gc(verbose = FALSE))

## ---- 2. write 10x MTX + barcode lists ---------------------------
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

writeMM(counts, file.path(input_dir, "matrix.mtx"))
write.table(data.frame(rownames(counts), rownames(counts), "Gene Expression"),
            file.path(input_dir, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(colnames(counts), file.path(input_dir, "barcodes.tsv"))
for (f in c("matrix.mtx", "features.tsv", "barcodes.tsv"))
  system2("gzip", c("-f", shQuote(file.path(input_dir, f))))

writeLines(keep, file.path(input_dir, "keep.txt"))
writeLines(ref,  file.path(input_dir, "ref.txt"))

message("wrote ", input_dir)
message(sprintf("  matrix: %d genes x %d cells (keep %d + ref %d)",
                nrow(counts), ncol(counts), length(keep), length(ref)))
