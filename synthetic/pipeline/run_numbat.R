#!/usr/bin/env Rscript
# Run Numbat on Standardized_Inputs/Numbat (counts + annotations + alleles).
library(numbat)
library(readr)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: run_numbat.R <Input_Dir> <Output_Dir>")

input_dir <- args[1]
output_dir <- args[2]

# 1. Load Data
message("Loading inputs...")
counts <- readRDS(file.path(input_dir, "counts.rds"))
df_annot <- read_tsv(file.path(input_dir, "annotations.tsv"), show_col_types = FALSE)
allele_file <- file.path(input_dir, "alleles.tsv.gz")

# Sanity checks
message("Sanity checks...")
message(paste0("  counts dims: ", nrow(counts), " genes x ", ncol(counts), " cells"))
message(paste0("  annotations rows: ", nrow(df_annot)))
if (!"cell_id" %in% colnames(df_annot)) stop("annotations.tsv missing column: cell_id")
if (!"group" %in% colnames(df_annot)) stop("annotations.tsv missing column: group")

common_cells <- intersect(colnames(counts), df_annot$cell_id)
message(paste0("  shared cells (counts ∩ annotations): ", length(common_cells)))
if (length(common_cells) == 0) {
  stop("No overlap between counts columns and annotations cell_id. Check barcode format.")
}

normal_cells <- df_annot$cell_id[df_annot$group == "Normal"]
message(paste0("  normal cells in annotations: ", length(normal_cells)))
normal_shared <- intersect(colnames(counts), normal_cells)
message(paste0("  normal cells with counts: ", length(normal_shared)))
if (length(normal_shared) == 0) {
  stop("No Normal cells found after intersecting annotations with counts.")
}

# 2. Define Reference
# We filter the annotation to get the reference cells for Numbat baseline
ref_cells <- df_annot$cell_id[df_annot$group == "Normal"]

# 3. Run Numbat
# Numbat automatically handles the separation of reference and query internally
# based on the 'lambdas_ref' or 'ref_annot' logic. 
# Here we calculate the reference expression profile (lambdas_ref) from Normal cells.

# Calculate reference expression profile
message("Calculating reference profile...")
ref_counts <- counts[, colnames(counts) %in% ref_cells]
lambdas_ref <- rowMeans(ref_counts)

ref_annot_internal <- data.frame(
  cell = colnames(ref_counts),
  group = 'Normal'
)
ref_internal <- aggregate_counts(ref_counts, ref_annot_internal)

## allele file 
df_allele <- read_tsv(allele_file, show_col_types = FALSE)
# Run
out <- run_numbat(
  count_mat = counts,
  lambdas_ref = ref_internal,
  df_allele = df_allele,
  genome = "hg38",
  out_dir = output_dir,
  ncores = 4, # Adjust based on your server
  plot = TRUE,
  common_diploid = FALSE
)
saveRDS(out, file.path(output_dir, "numbat_output.rds"))