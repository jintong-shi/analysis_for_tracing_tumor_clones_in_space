#!/usr/bin/env Rscript
# novaseq_08_g3_benign_ref_infcnv.R
#
# Benign-reference inferCNV on IDC_Grade3 (NovaSeq X Illumina spatial), mirroring
# opensource/benchmark_spatialCNV/realdata/highres/atera/benign_ref_strictqc_infcnv.R.
#
# Cells come from novaseq_07_g3_select_regions.py: epithelial by clustering
# (major_celltype == "epi", frozen res 0.3) AND inside a histology polygon --
#   group "normal" = the 4 benign polygons (n = 476, all but one in cluster 10)
#   group "cancer" = the IDC polygon        (n = 519,528)
#
# inferCNV parameters identical to the atera script:
#   cutoff = 0.1, denoise = TRUE, HMM = FALSE, analysis_mode = "cells",
#   cluster_by_groups = FALSE, write_phylo = TRUE, num_threads = 4,
#   chr_exclude = "chrM", gene_order_file = gencode_v38.
#
# Two deliberate departures, neither an analysis parameter:
#   - save_rds = FALSE. /data is 98% full (682G free); at 520k cells the per-step
#     RDS checkpoints would fill the disk. Only the final object is written.
#   - the heatmaps live in novaseq_09_g3_infcnv_heatmap.R: no conda env has both
#     infercnv and ComplexHeatmap, and the all-cell panel is capped at 50,000 cells
#     instead of NULL -- clustering 520k rows in ComplexHeatmap is not feasible.
#
# Env: stcnv-r-shared.
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_08_g3_benign_ref_infcnv.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
})

set.seed(42)

## ---- config -----------------------------------------------------
sample        <- "IDC_Grade3"
proj          <- "<DATA_ROOT>/high_res_ST/novaseq"
counts_h5     <- file.path(proj, "processed_data", "cells", paste0(sample, "_cells_filtered.h5"))
groups_csv    <- file.path(proj, "processed_data", paste0("infercnv_", sample),
                           paste0(sample, "_epi_groups.csv"))
output_dir    <- file.path(proj, "processed_data", paste0("infercnv_", sample))
gene_pos_file <- "<REF_ROOT>/gencode_v38_gene_pos.txt"

cutoff        <- 0.1
ncores        <- 4

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. counts + group labels -----------------------------------
message("Reading counts: ", counts_h5)
counts <- Read10X_h5(counts_h5)

groups <- read.csv(groups_csv, row.names = 1)
missing <- setdiff(rownames(groups), colnames(counts))
if (length(missing) > 0) {
  stop(length(missing), " selected cells are absent from the counts matrix.")
}
counts <- counts[, rownames(groups)]
message(sprintf("selected %d cells x %d genes", ncol(counts), nrow(counts)))
print(table(groups$group))

group.dat <- data.frame(row.names = rownames(groups), group = groups$group)

gene_pos <- read.delim(gene_pos_file, header = FALSE, stringsAsFactors = FALSE)
message(sprintf("genes shared with the gene-order file: %d / %d",
                length(intersect(rownames(counts), gene_pos[[1]])), nrow(counts)))

## ---- 2. benign-reference inferCNV -------------------------------
infcnv_dir <- file.path(output_dir, "infcnv")
dir.create(infcnv_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(infcnv_dir, "run.final.infercnv_obj"))) {
  message("Running benign-reference inferCNV (ref = normal) ...")
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    gene_order_file   = gene_pos_file,
    annotations_file  = group.dat,
    delim             = "\t",
    ref_group_names   = "normal",
    chr_exclude       = "chrM"
  )
  rm(counts); gc()
  infercnv::run(
    obj, cutoff = cutoff, out_dir = infcnv_dir,
    denoise = TRUE, HMM = FALSE,
    cluster_by_groups = FALSE,
    write_phylo = TRUE,
    analysis_mode = "cells",
    num_threads = ncores,
    no_plot = TRUE,
    save_rds = TRUE
  )
}
final_obj <- file.path(infcnv_dir, "run.final.infercnv_obj")
message("inferCNV done: ", final_obj)

message("Done. Plot with novaseq_09_g3_infcnv_heatmap.R (env stcnv-r-tuscan).")
