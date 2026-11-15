#!/usr/bin/env Rscript
# Run CONICSmat on Standardized_Inputs/CONICSmat (counts.rds + references.rds).
suppressPackageStartupMessages({
  library(CONICSmat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(zoo)
  library(pdist)
  library(data.table)
})

# ==========================================
# 1. Parse Arguments
# ==========================================
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Usage: run_conicsmat.R <Input_Dir> <Output_Dir> <Chrom_Arm_File> [Gene_Pos_File]")
}

input_dir      <- args[1]
output_dir     <- args[2]
chrom_arm_file <- args[3]
gene_pos_file  <- if (length(args) >= 4) args[4] else ""

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# ==========================================
# 2. Load Data
# ==========================================
message("-> Loading inputs...")
counts_file <- file.path(input_dir, "counts.rds")
ref_file    <- file.path(input_dir, "references.rds")

if (!file.exists(counts_file)) stop("Counts file missing.")
if (!file.exists(ref_file))    stop("Reference file missing.")
if (!file.exists(chrom_arm_file)) stop("Chrom arm file missing.")
if (nchar(gene_pos_file) > 0 && !file.exists(gene_pos_file)) {
  stop(paste("Gene position file missing:", gene_pos_file))
}

# Load Counts (Genes x Cells)
# CONICSmat expects standard matrix or dgCMatrix
counts <- readRDS(counts_file)
if (is(counts, "sparseMatrix")) counts <- as.matrix(counts)

# Load References
ref_barcodes <- readRDS(ref_file)
if (is.list(ref_barcodes)) ref_barcodes <- unlist(ref_barcodes)

# Load Regions
regions <- read.table(chrom_arm_file, sep="\t", row.names=1, header=TRUE)

message(sprintf("   Matrix: %d genes x %d spots", nrow(counts), ncol(counts)))
message(sprintf("   Reference Cells: %d", length(ref_barcodes)))

# ==========================================
# 3. Preprocessing (CPM -> Log2)
# ==========================================
message("-> Preprocessing...")

# CPM and Log2 transform (suva_expr)
total_counts <- colSums(counts)
cpm_matrix   <- sweep(counts, 2, total_counts, FUN="/") * 1e6
suva_expr    <- log2(cpm_matrix / 10 + 1)
suva_expr[is.na(suva_expr)] <- 0

# Get Gene Positions
# Using the archive version from your example to ensure stability
load_gene_pos <- function(fp) {
  gp <- data.table::fread(fp, data.table = FALSE)
  if (ncol(gp) < 4) stop("Gene position file must have at least 4 columns.")
  nm <- tolower(names(gp))
  pick <- function(cands) {
    idx <- which(nm %in% cands)
    if (length(idx) > 0) return(idx[1])
    return(NA_integer_)
  }
  gene_col <- pick(c("gene","hgnc_symbol","symbol"))
  chr_col  <- pick(c("chr","chromosome","chromosome_name"))
  start_col <- pick(c("start","start_position","txstart"))
  end_col   <- pick(c("end","end_position","txend"))
  if (any(is.na(c(gene_col, chr_col, start_col, end_col)))) {
    gene_col <- 1; chr_col <- 2; start_col <- 3; end_col <- 4
  }
  gene <- as.character(gp[[gene_col]])
  chr  <- as.character(gp[[chr_col]])
  start <- suppressWarnings(as.numeric(gp[[start_col]]))
  end   <- suppressWarnings(as.numeric(gp[[end_col]]))
  gene <- gsub('"', '', gene)
  chr  <- gsub('"', '', chr)
  chr  <- ifelse(grepl("^chr", chr, ignore.case = TRUE), chr, paste0("chr", chr))
  data.frame(
    chromosome_name = chr,
    hgnc_symbol = gene,
    start_position = start,
    end_position = end,
    stringsAsFactors = FALSE
  )
}

if (nchar(gene_pos_file) > 0) {
  message("-> Loading gene positions from local file...")
  gene_pos <- load_gene_pos(gene_pos_file)
} else {
  message("-> Fetching gene positions (Ensembl)...")
  gene_pos <- getGenePositions(rownames(suva_expr), ensembl_version="https://www.ensembl.org")
}

# Filter Matrix
suva_expr <- filterMatrix(suva_expr, gene_pos[,"hgnc_symbol"], minCells=5)

# ==========================================
# 4. Calculate CNV Matrix (Ref-based Manual Logic)
# ==========================================
message("-> Calculating CNV matrix (Reference subtraction)...")

expThresh <- 0.1
windowsize <- 101

# Identify Normal (Ref) vs Tumor
# Ensure ref_barcodes are actually in the matrix columns
ref_clean <- intersect(ref_barcodes, colnames(suva_expr))
tumor_clean <- setdiff(colnames(suva_expr), ref_clean)

if(length(ref_clean) < 2) stop("Not enough valid reference cells found in matrix.")

# Calculate Reference Vector
ref_vec <- rowMeans(suva_expr[, ref_clean])

# Subset Matrix to Ref + Tumor (ensure order)
gexp <- suva_expr[, c(ref_clean, tumor_clean)]

# Filter Low Expression Genes based on Reference and Overall
ref_vec <- ref_vec[which(ref_vec > expThresh)]
gexp    <- gexp[names(ref_vec), ]
gexp    <- gexp[which(rowMeans(gexp) > expThresh), ]

# Align ref_vec again
ref_vec <- ref_vec[rownames(gexp)]

# Clean Gene Positions for Ordering
gene_pos_clean <- gene_pos %>%
  mutate(
    chromosome_name = gsub("^chr", "", chromosome_name, ignore.case = TRUE),
    chr_i   = suppressWarnings(as.integer(as.character(chromosome_name))),
    start_i = map_dbl(start_position, ~ suppressWarnings(as.double(.x[1])))
  ) %>%
  filter(chr_i %in% 1:24, is.finite(start_i)) %>%
  arrange(chr_i, start_i) %>%
  select(-chr_i, -start_i)

colnames(gene_pos_clean)[2] <- 'gene'

# Align Matrix to Genomic Order
common_genes <- intersect(gene_pos_clean$gene, rownames(gexp))
gexp <- gexp[common_genes, ]
ref_vec <- ref_vec[common_genes]
gene_pos_final <- gene_pos_clean %>% filter(gene %in% common_genes)

# Calculate Ratio
rat <- gexp - ref_vec

# ==========================================
# 5. Smoothing (Rolling Average)
# ==========================================
message("-> Applying rolling average and clustering...")

# Rolling Mean
d <- apply(rat, 2, function(x) zoo::rollapply(x, mean, width = windowsize, align = "center"))

# Centering
d <- apply(d, 2, function(x) x - mean(x))

message(paste("   Distance matrix ready:", nrow(d), "genes x", ncol(d), "cells"))

# ==========================================
# 6. Clustering and Sorting
# ==========================================
# Indices in 'd'
idx_ref   <- 1:length(ref_clean)
idx_tumor <- (length(ref_clean) + 1):ncol(d)

# Cluster Tumor cells only
dmat <- as.dist(pdist(t(d[, idx_tumor])))
hc   <- hclust(dmat)

# Reorder: Ref first, then sorted Tumor
cellOrder <- c(idx_ref, idx_tumor[hc$order])
d_sorted  <- d[, cellOrder]

# ==========================================
# 7. Save Outputs
# ==========================================
message("-> Saving results...")

# Save the final matrix (Genes x Cells)
saveRDS(d_sorted, file.path(output_dir, "CONICS_track_matrix.rds"))

# Save the gene axis info (needed for plotting/evaluation alignment)
write.csv(gene_pos_final, file.path(output_dir, "CONICS_axis.csv"), row.names = FALSE)

message("-> Done.")
