#!/usr/bin/env Rscript
# visiumhd_20_cnvtools_prep.R
#
# Export the VisiumHD AP5 pseudo-spot (metacell) set as a 10x MTX directory so
# the shared pipeline scripts snake_pipeline/R/run_conicsmat.R and run_copykat.R
# can be run on it unmodified.
#
# Input is the r40_Smax20 pseudo-spot Seurat object (region-growing metacells,
# code/region_growing.R + processed_data/pseudo_spot/code/generate_pseudo_spot_seu.R),
# 18,060 genes x 4,721 pseudo-spots:
#   ref  = pure_type "Benign" (2,064 spots)
#   keep = pure_type "Tumor"  (2,158 spots)
# The 499 spots with status "Mixed" (pure_type NA) are dropped -- they are not
# assignable to either side, so they cannot serve as reference or as observation.
#
# Env: stcnv-r-shared.
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/visiumHD/code/visiumhd_20_cnvtools_prep.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

## ---- config -----------------------------------------------------
HD        <- "<DATA_ROOT>/high_res_ST/visiumHD"
seu_path  <- file.path(HD, "processed_data/pseudo_spot/epi_pseudo_r40_Smax20_seu.rds.gz")
input_dir <- file.path(HD, "processed_data/cnvtools/pseudo_r40_Smax20/input")

## ---- 1. counts + labels -----------------------------------------
message("reading pseudo-spot Seurat object (659 MB gz) ...")
seu <- readRDS(seu_path)

print(table(seu$pure_type, seu$status, useNA = "ifany"))

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
pure   <- seu$pure_type

keep <- colnames(seu)[!is.na(pure) & pure == "Tumor"]
ref  <- colnames(seu)[!is.na(pure) & pure == "Benign"]
if (length(ref) == 0L || length(keep) == 0L)
  stop("expected both Benign and Tumor pseudo-spots in pure_type")

message(sprintf("ref (Benign) = %d spots, keep (Tumor) = %d spots, dropped %d mixed",
                length(ref), length(keep), ncol(seu) - length(ref) - length(keep)))

counts <- counts[, c(keep, ref), drop = FALSE]
counts <- as(as(counts, "CsparseMatrix"), "generalMatrix")
rm(seu); invisible(gc(verbose = FALSE))

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
message(sprintf("  matrix: %d genes x %d spots (keep %d + ref %d)",
                nrow(counts), ncol(counts), length(keep), length(ref)))
