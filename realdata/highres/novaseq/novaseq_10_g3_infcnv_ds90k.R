#!/usr/bin/env Rscript
# novaseq_10_g3_infcnv_ds90k.R
#
# Same benign-reference inferCNV as novaseq_08_g3_benign_ref_infcnv.R, but with the
# cancer epithelial cells downsampled to 90,000. Written as a separate script so the
# all-cell version (08) stays untouched.
#
# Why downsample: inferCNV STEP 15 runs hclust(parallelDist(t(expr))) on the whole
# observation group. That distance matrix is n(n-1)/2 * 8 bytes -- 1.08 PB at 519,528
# cells, which is what killed the 08 run. At 90,000 cells it is 32 GB. For reference,
# the atera Xenium run (run4_benign_ref_strict_qc) did this with 84,247 cancer cells.
# hclust itself is not the limit: inferCNV resolves it to fastcluster::hclust, which
# has no 65,536-observation cap.
#
# Cells come from novaseq_07_g3_select_regions.py:
#   group "normal" = all 476 cells in the 4 benign polygons (reference, kept whole)
#   group "cancer" = 90,000 of the 519,528 cells in the IDC polygon, sampled within
#                    each frozen res-0.3 cluster in proportion to its size (seed 42)
#
# inferCNV parameters identical to the atera script and to 08:
#   cutoff = 0.1, denoise = TRUE, HMM = FALSE, analysis_mode = "cells",
#   cluster_by_groups = FALSE, write_phylo = TRUE, num_threads = 4,
#   chr_exclude = "chrM", gene_order_file = gencode_v38.
#   save_rds = FALSE because /data is 98% full.
#
# Env: stcnv-r-shared.
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_10_g3_infcnv_ds90k.R
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
output_dir    <- file.path(proj, "processed_data", paste0("infercnv_", sample, "_ds90k"))
gene_pos_file <- "<REF_ROOT>/gencode_v38_gene_pos.txt"

n_cancer      <- 90000
cutoff        <- 0.1
ncores        <- 4

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. downsample the cancer cells, stratified by cluster ------
groups <- read.csv(groups_csv, row.names = 1)
cancer <- groups[groups$group == "cancer", ]
normal <- groups[groups$group == "normal", ]

per_cluster <- table(cancer$cluster_frozen)
take <- floor(n_cancer * per_cluster / sum(per_cluster))
# hand out the rounding remainder to the largest clusters
short <- n_cancer - sum(take)
if (short > 0) {
  take[order(per_cluster, decreasing = TRUE)[seq_len(short)]] <-
    take[order(per_cluster, decreasing = TRUE)[seq_len(short)]] + 1
}

keep <- unlist(lapply(names(per_cluster), function(cl) {
  ids <- rownames(cancer)[cancer$cluster_frozen == cl]
  sample(ids, take[[cl]])
}), use.names = FALSE)

sel <- rbind(cancer[keep, ], normal)
message(sprintf("downsampled cancer %d -> %d; reference kept whole at %d",
                nrow(cancer), length(keep), nrow(normal)))
print(rbind(all = per_cluster, sampled = take))

write.csv(sel, file.path(output_dir, paste0(sample, "_epi_groups_ds90k.csv")))

## ---- 2. counts --------------------------------------------------
message("Reading counts: ", counts_h5)
counts <- Read10X_h5(counts_h5)

missing <- setdiff(rownames(sel), colnames(counts))
if (length(missing) > 0) {
  stop(length(missing), " selected cells are absent from the counts matrix.")
}
counts <- counts[, rownames(sel)]
message(sprintf("selected %d cells x %d genes", ncol(counts), nrow(counts)))
print(table(sel$group))

group.dat <- data.frame(row.names = rownames(sel), group = sel$group)

## ---- 3. benign-reference inferCNV -------------------------------
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
    subc
  )
}
message("inferCNV done: ", file.path(infcnv_dir, "run.final.infercnv_obj"))
message("Plot with: novaseq_09_g3_infcnv_heatmap.R ", infcnv_dir)
