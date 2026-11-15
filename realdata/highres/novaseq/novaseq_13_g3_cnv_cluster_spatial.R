#!/usr/bin/env Rscript
# Spatial map of the NovaSeq IDC_Grade3 ds90k cells, colored by CNV cluster.
#
# "CNV cluster" = the phase-1 dendrogram (ward.D2 on 50 PCs of the inferCNV matrix,
# all 90,000 cancer cells; built by altered_benign/phase1_precompute.R) cut at K.
# Clusters are renumbered in leaf order, so cluster 1..K read left-to-right in the
# app heatmap / top-to-bottom in the dendrogram.
#
# Layout follows atera_subclone_spatial_dimplot.R (x = col, y = row reversed).
#
# Run with tuscan R (ggplot2):
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_13_g3_cnv_cluster_spatial.R

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

## --labels <csv>: colour by an external per-cell label (columns barcode,label),
## e.g. the CNV k-means of novaseq_14_g3_cnv_kmeans.R, instead of the tree cut.
args    <- commandArgs(trailingOnly = TRUE)
lab_csv <- if ("--labels" %in% args) args[[which(args == "--labels") + 1L]] else NA_character_
tag     <- if (is.na(lab_csv)) "cnv_cluster" else   # e.g. "cnv_kmeans_k6"
             sub("^IDC_Grade3_ds90k_", "", sub("\\.csv$", "", basename(lab_csv)))
how     <- if (is.na(lab_csv)) "tree cut" else
           if (grepl("merged", tag)) "merged CNV clusters" else "k-means"

PROJ    <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir <- file.path(PROJ, "processed_data", "infercnv_IDC_Grade3_ds90k")
rds     <- file.path(run_dir, "app_outputs", "novaseq_idc_grade3_ds90k.rds")
grp_csv <- file.path(run_dir, "IDC_Grade3_epi_groups_ds90k.csv")
out_dir <- file.path(run_dir, "outputs")
out_png <- file.path(out_dir, sprintf("IDC_Grade3_ds90k_%s_spatial.png", tag))
out_csv <- file.path(out_dir, "IDC_Grade3_ds90k_cnv_cluster_k12.csv")  # tree cut only

K <- 12L

if (is.na(lab_csv)) {
  ## ---- cut the CNV tree, renumber in leaf order ----
  message("reading phase-1 RDS (~4.3 GB) ...")
  d      <- readRDS(rds)
  coords <- d$coords                     # barcode, row (y), col (x)
  cl     <- cutree(d$dendro, k = K)[d$leaf_bc]
  cl     <- factor(cl, levels = unique(cl), labels = seq_len(K))
  cnv    <- setNames(as.character(cl), d$leaf_bc)
} else {
  ## ---- external labels; coordinates from the groups table (no 4.3 GB read) ----
  g      <- read.csv(grp_csv, row.names = 1)
  coords <- data.frame(barcode = rownames(g), row = g$y, col = g$x,
                       stringsAsFactors = FALSE)
  lb     <- read.csv(lab_csv)            # barcode, label
  cnv    <- setNames(as.character(lb[[2]]), lb[[1]])
  K      <- length(unique(cnv))
  message(sprintf("labels from %s (k = %d)", lab_csv, K))
}

coords$cnv_cluster <- cnv[coords$barcode]
n_na <- sum(is.na(coords$cnv_cluster))
if (n_na > 0) message(sprintf("WARNING: %d cells have no label (dropped)", n_na))
coords <- coords[!is.na(coords$cnv_cluster), ]
coords$cnv_cluster <- factor(coords$cnv_cluster, levels = as.character(seq_len(K)))
message("CNV cluster counts:"); print(table(coords$cnv_cluster))

if (is.na(lab_csv)) {  # the k-means labels already live in their own CSV
  write.csv(coords[, c("barcode", "col", "row", "cnv_cluster")], out_csv, row.names = FALSE)
  on.exit(message("DONE ", out_csv), add = TRUE)
}

## ---- draw ----
pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4",
         "#F032E6", "#BFEF45", "#FABED4", "#469990", "#9A6324", "#000075")
pal        <- pal[seq_len(K)]          # K <= 12
names(pal) <- as.character(seq_len(K))
n_cl <- table(coords$cnv_cluster)
labs <- sprintf("%s  (%s)", names(n_cl), format(as.integer(n_cl), big.mark = ","))

p <- ggplot(coords, aes(x = col, y = row, color = cnv_cluster)) +
  geom_point(size = 0.35, stroke = 0) +
  scale_color_manual(values = pal, labels = labs, name = "CNV cluster", drop = FALSE) +
  scale_y_reverse() +
  coord_fixed() +
  labs(title = sprintf("NovaSeq IDC_Grade3 ds90k - CNV cluster spatial map, %s (k = %d, n = %d cells)",
                       how, K, nrow(coords))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

xr <- diff(range(coords$col)); yr <- diff(range(coords$row))
maxdim <- 11
if (xr >= yr) { w <- maxdim; h <- maxdim * yr / xr } else { h <- maxdim; w <- maxdim * xr / yr }
ggsave(out_png, p, width = w + 2.5, height = h, dpi = 300, limitsize = FALSE)
message("DONE ", out_png)

