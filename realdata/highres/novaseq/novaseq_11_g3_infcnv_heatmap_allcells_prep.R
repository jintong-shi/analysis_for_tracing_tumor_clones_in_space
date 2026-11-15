#!/usr/bin/env Rscript
# All-cell IDC_Grade3 ds90k inferCNV heatmap - prep (step 1/2).
#
# novaseq_09_g3_infcnv_heatmap.R samples 5,000 (or 50,000) cells and clusters the
# rows; the 50,000 version needs a 20 GB dist() and ~5-10 h, so it is unusable.
# This pair follows the atera all-cell pattern instead
# (high_res_ST/atera/code/atera_subclone_cnv_heatmap_allcells_{prep,plot}.R):
# no row clustering, rows grouped by an existing label, and rows bin-averaged so
# ComplexHeatmap's rasteriser has a few thousand rows instead of 90,000.
#
# Here the grouping label is cluster_frozen (res-0.3 expression clusters) from
# IDC_Grade3_epi_groups_ds90k.csv, i.e. 8 cancer blocks over 90,000 cells, plus the
# 476 benign reference cells as their own panel (kept whole, unbinned).
#
# Binning is done in this step, not at plot time as in atera, so the RDS handed to
# the render step is ~1 GB instead of ~14 GB (/data is 97% full).
# BIN = 12 (atera used 34): cluster_frozen is an expression clustering, not a
# CNV-derived subclone, so blocks are less homogeneous and a smaller bin keeps more
# within-block structure. 90,000 -> ~7,510 rows.
#
# Env: stcnv-r-shared (only env with infercnv).
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_11_g3_infcnv_heatmap_allcells_prep.R
# ------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(infercnv)
}))

## --labels <csv>: group rows by an external per-cell label (columns barcode,label)
## instead of cluster_frozen, e.g. the CNV k-means of novaseq_14_g3_cnv_kmeans.R.
args    <- commandArgs(trailingOnly = TRUE)
lab_csv <- if ("--labels" %in% args) args[[which(args == "--labels") + 1L]] else NA_character_
lb      <- if (is.na(lab_csv)) NULL else read.csv(lab_csv)      # barcode, label
tag     <- if ("--tag" %in% args) args[[which(args == "--tag") + 1L]] else
           if (is.na(lab_csv)) "allcells" else
             sprintf("allcells_kmeans_k%d", length(unique(lb[[2]])))

sample    <- "IDC_Grade3"
proj      <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir   <- file.path(proj, "processed_data", paste0("infercnv_", sample, "_ds90k"))
obj_path  <- file.path(run_dir, "infcnv", "run.final.infercnv_obj")
groups_csv<- file.path(run_dir, paste0(sample, "_epi_groups_ds90k.csv"))
out_dir   <- file.path(run_dir, "outputs")
out_rds   <- file.path(out_dir, paste0(sample, "_ds90k_", tag, "_heatmap_input.rds"))

BIN  <- 12L
SEED <- 42L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("loading inferCNV object (~1.3G) ...")
o     <- readRDS(obj_path)
expr  <- o@expr.data            # genes x cells
go    <- o@gene_order           # chr, start, stop
cells <- colnames(expr)

## ---- reference (benign) vs observation (cancer) ----
ref_idx <- unique(unlist(o@reference_grouped_cell_indices,   use.names = FALSE))
obs_idx <- unique(unlist(o@observation_grouped_cell_indices, use.names = FALSE))
if (length(ref_idx) == 0L || length(obs_idx) == 0L)
  stop("need both reference and observation cells in the inferCNV object")

## ---- group label per observation cell ----
if (is.na(lab_csv)) {
  grp <- read.csv(groups_csv, row.names = 1)
  cl  <- setNames(paste0("C", grp$cluster_frozen), rownames(grp))
} else {
  cl  <- setNames(paste0("C", lb[[2]]), lb[[1]])
  message("row grouping from ", lab_csv)
}
lab_obs  <- cl[cells[obs_idx]]
if (anyNA(lab_obs))
  stop(sprintf("%d observation cells have no group label", sum(is.na(lab_obs))))
message(sprintf("reference=%d  observation=%d", length(ref_idx), length(obs_idx)))
print(table(lab_obs))

## ---- genome-ordered, chromosome-split columns ----
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
chr_levels <- chr_levels[chr_levels %in% unique(as.character(go$chr))]
chr_genes  <- factor(as.character(go$chr), levels = chr_levels)
gene_ord   <- order(chr_genes, go$start)
chr_genes  <- chr_genes[gene_ord]

build <- function(idx) {
  mm <- t(expr[gene_ord, idx, drop = FALSE])   # cells x genes
  rownames(mm) <- cells[idx]
  mm
}
mat_ref <- build(ref_idx)
mat_obs <- build(obs_idx)
rm(expr, o); gc()

## ---- random row order within panel, then bin-average within each label ----
## Random first, so a bin is ~BIN random cells of one cluster (an unbiased mean of
## that cluster), never the BIN most similar cells - the latter would manufacture
## structure that is not in the data.
set.seed(SEED)
mat_ref <- mat_ref[sample(nrow(mat_ref)), , drop = FALSE]
ord     <- sample(nrow(mat_obs))
mat_obs <- mat_obs[ord, , drop = FALSE]
lab_obs <- lab_obs[ord]

lv    <- paste0("C", sort(unique(as.integer(sub("^C", "", lab_obs)))))
parts <- lapply(lv, function(l) {
  mm <- mat_obs[lab_obs == l, , drop = FALSE]
  g  <- ceiling(seq_len(nrow(mm)) / BIN)
  rowsum(mm, g) / as.numeric(table(g))          # per-bin row means
})
n_per_lab <- vapply(lv, function(l) sum(lab_obs == l), integer(1))
bin_lab   <- unlist(lapply(lv, function(l) rep(l, ceiling(n_per_lab[[l]] / BIN))))
mat_bin   <- do.call(rbind, parts)
message(sprintf("binned observation rows: %d -> %d (bin=%d)",
                nrow(mat_obs), nrow(mat_bin), BIN))

saveRDS(list(
  mat_ref   = mat_ref,
  mat_bin   = mat_bin,
  bin_lab   = bin_lab,
  lab_levels= lv,
  n_per_lab = n_per_lab,
  chr_genes = chr_genes,
  n_ref     = nrow(mat_ref),
  n_obs     = nrow(mat_obs),
  bin       = BIN
), out_rds)
message("DONE prep -> ", out_rds)
