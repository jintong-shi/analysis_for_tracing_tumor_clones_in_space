#!/usr/bin/env Rscript
# Run SCEVAN on Standardized_Inputs/SCEVAN (counts.rds + optional references.rds).
suppressPackageStartupMessages({
  library(SCEVAN)
  library(Matrix)
  library(dplyr)
  library(BiocParallel)
})

# ==========================================
# 1. Parse Arguments
# ==========================================
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: run_scevan.R <Input_Dir> <Output_Dir> [Threads]")
}

input_dir  <- args[1]
output_dir <- args[2]
# Default to 1 threads if not provided
n_threads  <- if (length(args) >= 3) as.numeric(args[3]) else 1

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# IMPORTANT: SCEVAN outputs intermediate files (plots, etc.) to the current directory.
setwd(output_dir)

# ==========================================
# 2. Load Data
# ==========================================
message("-> Loading inputs...")
counts_file <- file.path(input_dir, "counts.rds")
ref_file    <- file.path(input_dir, "references.rds")

if (!file.exists(counts_file)) stop(paste("Counts file missing:", counts_file))

# Load Count Matrix (Genes x Cells)
raw_counts <- readRDS(counts_file)

# Load References if available
norm_cells <- NULL
if (file.exists(ref_file)) {
  norm_cells <- readRDS(ref_file)
  if (is.list(norm_cells)) norm_cells <- unlist(norm_cells)
  message(sprintf("   Reference mode enabled: %d normal cells found.", length(norm_cells)))
} else {
  message("   Reference file not found. Running in reference-free mode.")
}

# ==========================================
# 3. Run SCEVAN
# ==========================================
message("-> Running SCEVAN pipelineCNA...")

run_serial <- function() {
  BiocParallel::register(SerialParam(), default = TRUE)
}

run_pipeline <- function(use_subclones = TRUE, force_serial = FALSE) {
  if (force_serial) {
    run_serial()
  }
  args <- list(
    count_mtx = raw_counts,
    par_cores = if (force_serial) 1 else n_threads,
    plotTree = TRUE,
    organism = "human",
    SUBCLONES = use_subclones
  )
  if (!is.null(norm_cells) && length(norm_cells) > 0) {
    args$norm_cell <- norm_cells
  }
  do.call(pipelineCNA, args)
}

primary_error <- NULL
results <- tryCatch(
  run_pipeline(use_subclones = TRUE),
  error = function(e) {
    primary_error <<- e
    NULL
  }
)

if (is.null(results)) {
  err_msg <- conditionMessage(primary_error)
  tsne_subclone_err <- grepl(
    "Error in data.frame\\(x = tsne\\$Y\\[, 1\\], y = tsne\\$Y\\[, 2\\], Subclones = pred\\)",
    err_msg
  ) || grepl("arguments imply differing number of rows: 0", err_msg, fixed = TRUE)
  open_port_err <- grepl("cannot find an open port", err_msg, fixed = TRUE)

  if (tsne_subclone_err) {
    message("-> Detected SCEVAN subclone/tSNE plotting error. Retrying with SUBCLONES=FALSE ...")
    results <- run_pipeline(use_subclones = FALSE)
    message("-> Retry completed with SUBCLONES=FALSE.")
  } else if (open_port_err) {
    message("-> Detected BiocParallel port allocation error. Retrying in serial mode ...")
    results <- run_pipeline(use_subclones = TRUE, force_serial = TRUE)
    message("-> Retry completed in serial mode.")
  } else {
    stop(primary_error)
  }
}

# ==========================================
# 4. Save Results
# ==========================================
message("-> Saving final RDS object...")
saveRDS(results, file.path(output_dir, "scevan_results.rds"))

message("-> Done. Output in: ", output_dir)
