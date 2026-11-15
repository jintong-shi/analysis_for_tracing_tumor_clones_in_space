#!/usr/bin/env Rscript
# novaseq_20_g3_cnvtools_prep.R
#
# Export the IDC_Grade3 ds90k cell set as a 10x MTX directory so the shared
# pipeline scripts snake_pipeline/R/run_conicsmat.R and run_copykat.R can be run
# on it unmodified.
#
# Cell set is the frozen ds90k downsample already used for inferCNV
# (novaseq_10_g3_infcnv_ds90k.R), so all three tools see the same cells:
#   ref  = group "normal" (476 cells, the 4 benign polygons)
#   keep = group "cancer" (90,000 cells, downsampled from the IDC polygon)
#
# Only IDC_Grade3 is in scope; DCIS_IDC_Grade1 and IDC_Grade2 are not run.
#
# Numbat is not included: the NovaSeq X download ships no BAM/CRAM (only
# Raw_Matrix_Files, segmentation contours and the OME-TIFF), and numbat needs a
# pileup over aligned reads for allele counts.
#
# Env: stcnv-r-shared.
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_20_g3_cnvtools_prep.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

## ---- config -----------------------------------------------------
sample     <- "IDC_Grade3"
proj       <- "<DATA_ROOT>/high_res_ST/novaseq"
counts_h5  <- file.path(proj, "processed_data", "cells", paste0(sample, "_cells_filtered.h5"))
groups_csv <- file.path(proj, "processed_data", paste0("infercnv_", sample, "_ds90k"),
                        paste0(sample, "_epi_groups_ds90k.csv"))
input_dir  <- file.path(proj, "processed_data", "cnvtools",
                        paste0(sample, "_ds90k"), "input")

## ---- 1. counts + group labels -----------------------------------
message("reading counts (3.1 GB h5): ", counts_h5)
counts <- Read10X_h5(counts_h5)

groups  <- read.csv(groups_csv, row.names = 1)
missing <- setdiff(rownames(groups), colnames(counts))
if (length(missing) > 0)
  stop(length(missing), " selected cells are absent from the counts matrix.")

print(table(groups$group))

keep <- rownames(groups)[groups$group == "cancer"]
ref  <- rownames(groups)[groups$group == "normal"]
if (length(ref) == 0L || length(keep) == 0L)
  stop("expected both a 'cancer' and a 'normal' group in ", groups_csv)

counts <- counts[, c(keep, ref), drop = FALSE]
counts <- as(as(counts, "CsparseMatrix"), "generalMatrix")
invisible(gc(verbose = FALSE))

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
