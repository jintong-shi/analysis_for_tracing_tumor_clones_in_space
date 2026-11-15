#!/usr/bin/env Rscript
# run_infercnv.R stopped at infercnv::run(up_to_step = 17): produces the HMM
# pred_cnv_genes.dat without the slow HMM tail.
library(infercnv)
library(Matrix)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: run_infercnv_hmm17.R <Input_Dir> <Output_Dir> <Gene_Order_File> [--HMM TRUE|FALSE]")
}

input_dir  <- args[1]
output_dir <- args[2]
gene_order_file <- args[3]
hmm_enabled <- TRUE

# Optional flag: --HMM TRUE|FALSE
if (length(args) > 3) {
  hmm_idx <- which(args == "--HMM")
  if (length(hmm_idx) > 0) {
    if (hmm_idx[1] == length(args)) stop("Missing value after --HMM (use TRUE or FALSE).")
    hmm_val <- tolower(args[hmm_idx[1] + 1])
    if (hmm_val %in% c("1", "true", "t", "yes", "y")) {
      hmm_enabled <- TRUE
    } else if (hmm_val %in% c("0", "false", "f", "no", "n")) {
      hmm_enabled <- FALSE
    } else {
      stop(paste0("Invalid --HMM value: ", args[hmm_idx[1] + 1], ". Use TRUE or FALSE."))
    }
  }
}

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==========================================
# 2. Load Data
# ==========================================
message("Loading counts and annotations...")
counts_file <- file.path(input_dir, "counts.rds")
annot_file  <- file.path(input_dir, "annotations.txt")

if (!file.exists(counts_file)) stop(paste("Counts file missing:", counts_file))
if (!file.exists(annot_file)) stop(paste("Annotation file missing:", annot_file))

counts <- readRDS(counts_file)
annot <- read.table(annot_file)

# ==========================================
# 3. Create InferCNV Object
# ==========================================
# Use the formatted TSV file here
obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annot,
  gene_order_file   = gene_order_file,
  delim             = "\t",
  ref_group_names = c("Normal")
)

# ==========================================
# 4. Run InferCNV (up to step 19 -> produces HMM pred_cnv_genes.dat)
# ==========================================
message(paste0("Running inferCNV with HMM = ", hmm_enabled, " (up_to_step = 17)"))
res <- infercnv::run(
  obj,
  cutoff = 0.1,
  out_dir = output_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = hmm_enabled,
  num_threads = 4,
  HMM_report_by = 'cell',
  up_to_step = 17
)

file.create(file.path(output_dir, "infercnv.OK"))
quit(save = "no", status = 0)