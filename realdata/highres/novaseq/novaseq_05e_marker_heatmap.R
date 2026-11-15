#!/usr/bin/env Rscript
# Cluster x top-marker heatmap for res0.3 -- the evidence behind the compartment call.
#
# Reads the capped matrix novaseq_05d_cluster_markers.R step 'extract' already wrote
# (2,000 cells per cluster), so the 7 GB object is not opened again, and the top-10
# table step 'test' wrote. Each gene's per-cluster mean is z-scored across clusters,
# which is the only way twelve clusters fit one colour scale when CPB1 and CD3D
# differ by orders of magnitude.
#
# Genes are ordered by the cluster they mark, so the block diagonal is there by
# construction and is NOT the evidence. What the figure actually shows is how sharp
# each block is and which blocks are shared: the row dendrogram (ward.D2 over all the
# marker genes) says which clusters are really the same thing.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05e_marker_heatmap.R

suppressWarnings(suppressMessages({
  library(data.table); library(Matrix); library(pheatmap)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
TOP_N  <- 10L

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

b   <- readRDS(file.path(OUT, sprintf("%s_res%s_markers_input.rds", SAMPLE, RES)))
top <- fread(file.path(OUT, sprintf("%s_res%s_top%d_markers.csv", SAMPLE, RES, TOP_N)))
ann <- fread(file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES)))
ann[, cluster := as.character(cluster)]
top[, cluster := as.character(cluster)]

cl_lv <- as.character(sort(as.integer(ann$cluster)))
cl    <- factor(b$idents, levels = cl_lv)
log("cells: ", length(cl), "  clusters: ", nlevels(cl))

# genes in cluster order; a gene that marks several clusters is kept where it ranks best
top[, ord := as.integer(factor(cluster, levels = cl_lv)) * 1000L + rank]
genes <- top[order(ord)][!duplicated(gene)]
log("marker genes: ", nrow(genes), " unique of ", nrow(top), " cluster-gene pairs")

X <- b$data[genes$gene, , drop = FALSE]
M <- t(sapply(cl_lv, function(g) Matrix::rowMeans(X[, cl == g, drop = FALSE])))
Z <- t(scale(t(t(M))))                       # z-score each GENE across clusters
Z <- t(Z); Z[!is.finite(Z)] <- 0             # back to clusters x genes
lim <- quantile(abs(Z), 0.99, na.rm = TRUE)

fwrite(data.table(cluster = rownames(Z), as.data.table(round(Z, 4))),
       file.path(OUT, sprintf("%s_res%s_marker_heatmap_z.csv", SAMPLE, RES)))

lab <- sprintf("%s  %s%s  (n=%s)", ann$cluster[match(cl_lv, ann$cluster)],
               ann$major_celltype[match(cl_lv, ann$cluster)],
               ifelse(ann$ambiguous[match(cl_lv, ann$cluster)], "*", ""),
               formatC(ann$n_cells[match(cl_lv, ann$cluster)], big.mark = ",", format = "d"))
rownames(Z) <- lab

arow <- data.frame(call = factor(ann$major_celltype[match(cl_lv, ann$cluster)],
                                 levels = c("epi", "imm", "stroma")),
                   row.names = lab)
acol <- data.frame(marks = factor(genes$cluster, levels = cl_lv), row.names = genes$gene)
acols <- list(call  = c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868"),
              marks = setNames(rep(grDevices::hcl.colors(6, "Set2"), length.out = length(cl_lv)), cl_lv))

# pdf to read at size, png because st-a100 has no pdf rasteriser to preview with
for (ext in c("pdf", "png"))
pheatmap(Z, cluster_rows = TRUE, clustering_method = "ward.D2", cluster_cols = FALSE,
         gaps_col = cumsum(rle(as.character(genes$cluster))$lengths),
         annotation_row = arow, annotation_col = acol, annotation_colors = acols,
         breaks = seq(-lim, lim, length.out = 101),
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         fontsize_row = 7, fontsize_col = 5, angle_col = 90, border_color = NA,
         main = sprintf("%s res %s -- top %d markers per cluster, per-cluster mean z-scored across clusters\n* = compartment call flagged ambiguous; column blocks are the cluster each gene marks",
                        SAMPLE, RES, TOP_N),
         filename = file.path(OUT, sprintf("%s_res%s_marker_heatmap.%s", SAMPLE, RES, ext)),
         width = 16, height = 4.6)
log("done -> ", OUT)
