#!/usr/bin/env Rscript
# Subclone CNV prep (step 1/2) - extract from an inferCNV object.
#
# Cuts the stored ward.D2 hclust over the tumor (observation) spots into k
# subclones (cutree), then saves a compact RDS holding everything the plot step
# needs: the genome-ordered obs CNV matrix, per-spot subclone labels, the
# per-subclone pseudobulk profile (+ a reference pseudobulk row), and the
# chromosome column split. Mirrors atera_subclone_cnv_heatmap_allcells_prep.R.
#
# Subclones are numbered by their left-to-right order along the dendrogram, so
# subclone_1/2/3 read in the same order as the heatmap dendrogram leaves.
#
# Run with shared R (has infercnv):
#   <CONDA_ROOT>/miniconda3/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/visiumHD/code/subclone_cnv_prep.R <obj_path> <k> <out_rds> <label>
#
# Args:
#   obj_path : path to run.final.infercnv_obj
#   k        : number of subclones (cutree k)
#   out_rds  : output RDS path
#   label    : dataset label used in plot titles

suppressWarnings(suppressMessages({
  library(infercnv)
}))

args     <- commandArgs(trailingOnly = TRUE)
obj_path <- args[1]
k        <- as.integer(args[2])
out_rds  <- args[3]
label    <- args[4]

o    <- readRDS(obj_path)
expr <- o@expr.data                     # genes x cells (residual CNV expression)
go   <- o@gene_order                    # chr, start, stop
cells<- colnames(expr)

obs_idx <- o@observation_grouped_cell_indices[[1]]
ref_idx <- o@reference_grouped_cell_indices[[1]]
obs_cells <- cells[obs_idx]
ref_cells <- cells[ref_idx]

## ---- subclones from the stored ward.D2 hclust over the observation spots ----
hc <- o@tumor_subclusters$hc[["all_observations"]]
stopifnot(setequal(hc$labels, obs_cells))
cl <- cutree(hc, k = k)                 # named by hc$labels, ids 1..k (arbitrary)

# relabel clusters by first appearance along the dendrogram (left -> right)
leaf_order   <- hc$labels[hc$order]
first_seen   <- unique(cl[leaf_order])  # cluster ids in dendrogram order
relabel      <- setNames(seq_along(first_seen), first_seen)
sub_num      <- relabel[as.character(cl)]
names(sub_num) <- names(cl)
sub          <- setNames(paste0("subclone_", sub_num), names(sub_num))
sub          <- sub[obs_cells]          # align to obs_cells order

message(sprintf("[%s] k=%d subclone breakdown:", label, k))
print(table(sub))

## ---- genome-ordered, chromosome-split columns ----
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
chr_levels <- chr_levels[chr_levels %in% unique(as.character(go$chr))]
chr_genes  <- factor(as.character(go$chr), levels = chr_levels)
gene_ord   <- order(chr_genes, go$start)
chr_genes  <- chr_genes[gene_ord]

# obs matrix: spots x genes (genome-ordered)
mat <- t(expr[gene_ord, obs_cells, drop = FALSE])

## ---- pseudobulk: mean CNV profile per subclone (+ reference) ----
sub_levels <- paste0("subclone_", seq_len(k))
pb_sub <- t(sapply(sub_levels, function(s) {
  colMeans(mat[sub == s, , drop = FALSE])
}))
ref_pb <- colMeans(t(expr[gene_ord, ref_cells, drop = FALSE]))
ref_name <- names(o@reference_grouped_cell_indices)[1]   # "Benign" / "Normal"
pb <- rbind(pb_sub, ref_pb)
rownames(pb) <- c(sub_levels, ref_name)

saveRDS(list(
  mat        = mat,          # obs spots x genes
  sub        = sub,          # named per-spot subclone label (obs_cells order)
  pb         = pb,           # (k subclones + ref) x genes pseudobulk
  ref_name   = ref_name,
  chr_genes  = chr_genes,    # per-column chromosome factor
  obs_cells  = obs_cells,
  k          = k,
  label      = label
), out_rds)
message("DONE prep -> ", out_rds,
        sprintf("  (obs=%d x genes=%d, k=%d)", nrow(mat), ncol(mat), k))
