#!/usr/bin/env Rscript
# One-time builder for input_data/strataMap/ (supp_strataMap.R).
#
# StrataMap (the Illumina platform behind the "novaseq" folder names), breast
# IDC grade 3. Three extracts, all small enough to sit in git-adjacent storage:
#
#   strataMap_supp_cells.csv     715,479 cells x 9 columns, pulled straight out
#     of the res0.3 annotation table: x, y, umap1, umap2, nCount, nFeature,
#     area_um2, cluster, celltype. Everything the spatial/UMAP panels plot.
#   strataMap_marker_avgexpr.csv    the res0.3 rows of the marker average-expression
#     table (16 markers x 12 clusters), for the EPCAM-by-cluster dot plot.
#   strataMap_cluster_scores.csv    the res0.3 annotation table: the epithelial /
#     immune / stromal module scores the epi-vs-rest call was actually made from,
#     plus the winner, its margin over the runner-up, and the ambiguous flag.
#   strataMap_percell_cnv.csv    the per-cell CNV panel, already reduced to what
#     is drawn: 1500 display cells x 281 10Mb bins for each of the three tools,
#     z-scored across bins, in the shared hclust cell order. Recomputed here the
#     same way high_res_ST/code/compare_cnvtools_similarity.R does it (same seed,
#     same N_SHOW, same ward.D2 on the inferCNV profiles), so the panel reproduces
#     results/cnvtool_similarity/percell_heatmap_novaseq_ds20k_10mb.pdf rather than
#     copying it. Loading the three tools' matrices takes a few minutes and ~30 GB
#     of reads, which is exactly why the result is frozen here.
#
# Not a figure script; kept for provenance. Run once, from anywhere:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_strataMap_extracts.R

suppressWarnings(suppressMessages(library(data.table)))

PROJ  <- "<DATA_ROOT>"
CLUST <- file.path(PROJ, "high_res_ST/novaseq/processed_data/clustering/IDC_Grade3")
H5    <- file.path(PROJ, "high_res_ST/novaseq/processed_data/cells/IDC_Grade3_cells_filtered.h5")
PLATFORM <- "novaseq_ds20k"; BIN_MB <- 10
N_SHOW   <- 1500L; SEED <- 1234L          # must match compare_cnvtools_similarity.R
TOOLS    <- c("inferCNV", "CopyKAT", "CONICSmat")

args <- commandArgs(trailingOnly = FALSE)
SELF <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
OUT  <- file.path(dirname(dirname(SELF)), "strataMap")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. per-cell table -------------------------------------------------------
a <- fread(file.path(CLUST, "IDC_Grade3_annotated_cells_res0.3.csv"))
setnames(a, "V1", "cell")
cells <- a[, .(cell, x, y, umap1, umap2, nCount = nCount_spatial,
               nFeature = nFeature_spatial, area_um2,
               cluster = as.character(cluster_frozen), celltype = major_celltype)]

## ---- 1b. per-cell EPCAM, pulled out of the raw count matrix ------------------
# One gene out of a 37,372 x 715,479 CSC matrix with 1.23 billion non-zeros. The
# indices array is scanned in chunks to find that gene's row, which costs seconds;
# the alternative source is the 7 GB annotated Seurat object. Counts become
# Seurat's LogNormalize (log1p of counts-per-10k) using the nCount already in the
# table, so the values match what a Seurat FeaturePlot would show.
suppressWarnings(suppressMessages(library(hdf5r)))
GENE <- "EPCAM"; CHUNK <- 5e7
h5   <- H5File$new(H5, "r")
row0 <- which(h5[["matrix/features/name"]]$read() == GENE) - 1L   # 0-based
if (length(row0) != 1L) stop(GENE, " is not a unique row of ", H5)
bc   <- h5[["matrix/barcodes"]]$read()
ptr  <- h5[["matrix/indptr"]]$read()                              # 0-based, length ncol+1
di   <- h5[["matrix/indices"]]; dd <- h5[["matrix/data"]]
nnz  <- di$dims
hit_at <- integer(0); hit_val <- numeric(0)
for (st in seq(1, nnz, by = CHUNK)) {
  en <- min(st + CHUNK - 1, nnz)
  ii <- di$read(args = list(st:en))
  k  <- which(ii == row0)
  if (length(k)) {
    hit_at  <- c(hit_at,  as.integer(k + st - 1L))
    hit_val <- c(hit_val, dd$read(args = list((st:en)[k])))
  }
}
h5$close_all()
cnt <- setNames(numeric(length(bc)), bc)
cnt[findInterval(hit_at - 1L, ptr)] <- hit_val
cells[, EPCAM := { n <- cnt[cell]; log1p(n / nCount * 1e4) }]
message(sprintf("%s: %d/%d cells non-zero, max lognorm %.2f",
                GENE, sum(cells$EPCAM > 0), nrow(cells), max(cells$EPCAM)))

fwrite(cells, file.path(OUT, "strataMap_supp_cells.csv"))
message(sprintf("cells: %d, clusters %s, celltypes %s", nrow(cells),
                paste(sort(unique(cells$cluster)), collapse = ","),
                paste(sort(unique(cells$celltype)), collapse = ",")))

## ---- 2. marker averages, res 0.3 only ---------------------------------------
m <- fread(file.path(CLUST, "IDC_Grade3_marker_avgexpr.csv"))
m <- m[resolution == 0.3]
fwrite(m, file.path(OUT, "strataMap_marker_avgexpr.csv"))
message(sprintf("markers: %d clusters x %d genes", nrow(m), ncol(m) - 2L))

## ---- 2b. the module scores the compartment call was made from ---------------
an <- fread(file.path(CLUST, "IDC_Grade3_annotation_res0.3.csv"))
an[, cluster := as.character(cluster)]
fwrite(an, file.path(OUT, "strataMap_cluster_scores.csv"))
message(sprintf("cluster scores: %d clusters, %d flagged ambiguous (%.0f%% of cells)",
                nrow(an), sum(an$ambiguous), 100 * sum(an$pct_cells[an$ambiguous]) /
                  sum(an$pct_cells)))

## ---- 3. per-cell CNV, the three tools on one shared cell order ---------------
source(file.path(PROJ, "high_res_ST/code/cnv_bin_utils.R"))
M     <- load_tool_bins(PLATFORM, BIN_MB)
keep  <- Reduce(intersect, c(list(M$keep), lapply(M[TOOLS], colnames)))
bins  <- Reduce(intersect, lapply(M[TOOLS], rownames))
bins  <- bins[order(as.integer(sub(":.*", "", bins)), as.integer(sub(".*:", "", bins)))]
mats  <- lapply(M[TOOLS], function(x) x[bins, keep]); rm(M); gc(FALSE)
message(sprintf("cnv: %d cells x %d bins", length(keep), length(bins)))

set.seed(SEED)
show <- if (length(keep) > N_SHOW) sort(sample(seq_along(keep), N_SHOW)) else seq_along(keep)
show <- show[hclust(dist(t(mats$inferCNV[, show])), method = "ward.D2")$order]

SC <- do.call(rbind, lapply(TOOLS, function(t) t(mats[[t]][, show, drop = FALSE])))
Z  <- t(scale(t(SC))); Z[!is.finite(Z)] <- 0
out <- data.table(tool = rep(TOOLS, each = length(show)),
                  cell = rep(keep[show], length(TOOLS)),
                  row_in_tool = rep(seq_along(show), length(TOOLS)))
fwrite(cbind(out, as.data.table(round(Z, 4))), file.path(OUT, "strataMap_percell_cnv.csv"))
message(sprintf("percell cnv: %d rows (%d cells x %d tools) x %d bins",
                nrow(Z), length(show), length(TOOLS), ncol(Z)))
message("done -> ", OUT)
