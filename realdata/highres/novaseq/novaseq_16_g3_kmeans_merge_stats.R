#!/usr/bin/env Rscript
# Evidence for merging the k=12 CNV k-means clusters of the ds90k run.
#
# Two similarity matrices over the 12 clusters:
#   CNV      - Pearson correlation between mean inferCNV profiles (centred on 1),
#              i.e. do two clusters carry the same gains and losses.
#   spatial  - cosine similarity of their cell densities over a 200-unit tissue
#              grid (do they sit in the same place), plus the distance between
#              their spatial centroids as a share of tissue diagonal (are they at
#              least adjacent).
# A merge wants high CNV correlation AND either high spatial overlap or a small
# centroid distance; CNV-similar but far-apart clusters are left alone.
#
# Writes the cluster mean profiles so a merge can be applied without reloading
# the 4.2 GB matrix.
#
# Env: stcnv-r-shared.
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_16_g3_kmeans_merge_stats.R

suppressWarnings(suppressMessages({
  library(infercnv)
}))

sample   <- "IDC_Grade3"
PROJ     <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir  <- file.path(PROJ, "processed_data", paste0("infercnv_", sample, "_ds90k"))
obj_path <- file.path(run_dir, "infcnv", "run.final.infercnv_obj")
grp_csv  <- file.path(run_dir, paste0(sample, "_epi_groups_ds90k.csv"))
out_dir  <- file.path(run_dir, "outputs")
lab_csv  <- file.path(out_dir, paste0(sample, "_ds90k_cnv_kmeans_k12.csv"))
cent_csv <- file.path(out_dir, paste0(sample, "_ds90k_cnv_kmeans_k12_centroids.csv"))

GRID <- 200                                    # tissue grid step, coordinate units

lb  <- read.csv(lab_csv)
cl  <- setNames(as.integer(lb[[2]]), lb[[1]])
lv  <- sort(unique(cl))

## ---- CNV: mean profile per cluster ----
message("loading inferCNV object ...")
o        <- readRDS(obj_path)
obs_idx  <- unique(unlist(o@observation_grouped_cell_indices, use.names = FALSE))
cells    <- colnames(o@expr.data)[obs_idx]
go       <- o@gene_order
chr_lv   <- paste0("chr", c(1:22, "X", "Y"))
chr_lv   <- chr_lv[chr_lv %in% unique(as.character(go$chr))]
gene_ord <- order(factor(as.character(go$chr), levels = chr_lv), go$start)
e <- o@expr.data[gene_ord, obs_idx, drop = FALSE]
rm(o); invisible(gc())

lab  <- cl[cells]
cent <- vapply(lv, function(l) rowMeans(e[, lab == l, drop = FALSE]),
               numeric(nrow(e)))
colnames(cent) <- paste0("C", lv)
rm(e); invisible(gc())
write.csv(cent, cent_csv)

cnv_cor <- cor(cent - 1)                       # centred: correlation of the deviations

## ---- spatial: density over a tissue grid + centroid distance ----
g  <- read.csv(grp_csv, row.names = 1)
g  <- g[names(cl), ]
bx <- paste(floor(g$x / GRID), floor(g$y / GRID))
tab <- table(bx, cl)                           # bins x clusters
tab <- sweep(tab, 2, colSums(tab), "/")        # each cluster -> density over bins
cs  <- crossprod(tab) / outer(sqrt(colSums(tab^2)), sqrt(colSums(tab^2)))
dimnames(cs) <- list(paste0("C", lv), paste0("C", lv))

ctr  <- t(vapply(lv, function(l) c(mean(g$x[cl == l]), mean(g$y[cl == l])), numeric(2)))
diag_len <- sqrt(diff(range(g$x))^2 + diff(range(g$y))^2)
ctr_d <- as.matrix(dist(ctr)) / diag_len
dimnames(ctr_d) <- dimnames(cs)

message("\n== CNV correlation between cluster mean profiles ==")
print(round(cnv_cor, 2))
message("\n== spatial overlap (cosine of grid densities) ==")
print(round(cs, 2))
message("\n== spatial centroid distance (share of tissue diagonal) ==")
print(round(ctr_d, 2))

message("\n== hclust of clusters on CNV profile (1 - correlation) ==")
hc <- hclust(as.dist(1 - cnv_cor), method = "average")
print(hc$labels[hc$order])
for (k in 2:8)
  message(sprintf("  cut k=%d: %s", k,
                  paste(sprintf("%s->%d", hc$labels, cutree(hc, k)), collapse = " ")))
message("\nwrote ", cent_csv)
