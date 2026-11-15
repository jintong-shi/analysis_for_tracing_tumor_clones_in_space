#!/usr/bin/env Rscript
# k-means on the inferCNV matrix of the NovaSeq IDC_Grade3 ds90k run.
#
# Groups cells by their CNV profile alone: the input is the inferCNV residual
# expression matrix (6,318 genes x 90,000 observation cells) with no dimension
# reduction and no expression-clustering label anywhere in the pipeline. The
# res-0.3 cluster_frozen labels are NOT used -- they are an expression clustering
# and were only ever the heatmap's row grouping.
#
# Reference (benign) cells are excluded: they are the baseline the matrix is
# already expressed relative to, not observations to be clustered.
#
# Env: stcnv-r-shared (needs infercnv for the S4 class definition).
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_14_g3_cnv_kmeans.R

suppressWarnings(suppressMessages({
  library(infercnv)
}))

args <- commandArgs(trailingOnly = TRUE)
K    <- if ("--k" %in% args) as.integer(args[[which(args == "--k") + 1L]]) else 12L

sample   <- "IDC_Grade3"
PROJ     <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir  <- file.path(PROJ, "processed_data", paste0("infercnv_", sample, "_ds90k"))
obj_path <- file.path(run_dir, "infcnv", "run.final.infercnv_obj")
out_dir  <- file.path(run_dir, "outputs")
out_csv  <- file.path(out_dir, sprintf("%s_ds90k_cnv_kmeans_k%d.csv", sample, K))

SEED  <- 42L
NSTART<- 2L
ITMAX <- 50L
## Provenance of the stored k=12 labels (2026-09-24). The CSV on disk was written
## by an earlier launch of this script with NSTART = 5, ITMAX = 100 (101 iterations,
## 11 h, still not converged); the constants above are the shorter settings used for
## everything downstream. So this file as committed does NOT byte-reproduce
## ..._cnv_kmeans_k12.csv. The two runs agree on the 3-way merge: every cluster of
## the long run is 98.6-100% pure for one merged M, 99.53% of cells keep their M,
## and the merge map of novaseq_17 is unchanged. Deliberately not resynced --
## novaseq_15..19 outputs are from the shorter run.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("loading inferCNV object (~1.3G) ...")
o       <- readRDS(obj_path)
obs_idx <- unique(unlist(o@observation_grouped_cell_indices, use.names = FALSE))
cells   <- colnames(o@expr.data)[obs_idx]

## genome order, so the matrix columns are comparable to the heatmap
go         <- o@gene_order
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
chr_levels <- chr_levels[chr_levels %in% unique(as.character(go$chr))]
gene_ord   <- order(factor(as.character(go$chr), levels = chr_levels), go$start)

## Drop the inferCNV object before transposing: it carries count.data and the HMM
## slots as well, and the node is memory-tight (jobs here have been OOM-killed).
e <- o@expr.data[gene_ord, obs_idx, drop = FALSE]
rm(o); invisible(gc())
m <- t(e); rownames(m) <- cells
rm(e); invisible(gc())
message(sprintf("matrix: %d cells x %d genes (%.1f GB)",
                nrow(m), ncol(m), as.numeric(object.size(m)) / 2^30))

set.seed(SEED)
t0 <- proc.time()[["elapsed"]]
km <- kmeans(m, centers = K, iter.max = ITMAX, nstart = NSTART, algorithm = "Lloyd")
message(sprintf("kmeans k=%d done in %.0f s (%d iterations, converged=%s)",
                K, proc.time()[["elapsed"]] - t0, km$iter, km$ifault == 0L))
message(sprintf("between_SS / total_SS = %.1f%%", 100 * km$betweenss / km$totss))
print(table(km$cluster))

write.csv(data.frame(barcode = names(km$cluster), kmeans = unname(km$cluster)),
          out_csv, row.names = FALSE)
message("DONE ", out_csv)
