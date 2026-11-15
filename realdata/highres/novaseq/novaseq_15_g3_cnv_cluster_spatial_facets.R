#!/usr/bin/env Rscript
# Per-cluster spatial facets for a CNV clustering of the NovaSeq IDC_Grade3 ds90k run.
#
# One small map per cluster: that cluster in its colour from the combined map
# (novaseq_13_g3_cnv_cluster_spatial.R), every other cell in light grey for tissue
# context. Reads the same labels CSV (barcode, label) and takes coordinates from the
# groups table, so it never touches the 4.3 GB phase-1 RDS.
#
# Run with tuscan R (ggplot2):
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_15_g3_cnv_cluster_spatial_facets.R \
#     --labels <outputs>/IDC_Grade3_ds90k_cnv_kmeans_k6.csv

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

args    <- commandArgs(trailingOnly = TRUE)
lab_csv <- args[[which(args == "--labels") + 1L]]
tag     <- sub("^IDC_Grade3_ds90k_", "", sub("\\.csv$", "", basename(lab_csv)))

PROJ    <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir <- file.path(PROJ, "processed_data", "infercnv_IDC_Grade3_ds90k")
grp_csv <- file.path(run_dir, "IDC_Grade3_epi_groups_ds90k.csv")
out_png <- file.path(run_dir, "outputs", sprintf("IDC_Grade3_ds90k_%s_spatial_facets.png", tag))

g      <- read.csv(grp_csv, row.names = 1)
lb     <- read.csv(lab_csv)
cnv    <- setNames(as.character(lb[[2]]), lb[[1]])
coords <- data.frame(barcode = rownames(g), row = g$y, col = g$x,
                     stringsAsFactors = FALSE)
coords$cnv_cluster <- cnv[coords$barcode]
coords <- coords[!is.na(coords$cnv_cluster), ]

n_cl <- table(coords$cnv_cluster)
n_cl <- n_cl[order(as.integer(names(n_cl)))]   # 1,2,..,12 -- not 1,10,11,12,2,..
K    <- length(n_cl)
lab  <- function(cl) sprintf("%s  (n = %s)", cl,
                             format(as.integer(n_cl[cl]), big.mark = ","))
coords$panel <- factor(lab(coords$cnv_cluster), levels = lab(names(n_cl)))
message(sprintf("%s: k = %d, %d cells", tag, K, nrow(coords)))

pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4",
         "#F032E6", "#BFEF45", "#FABED4", "#469990", "#9A6324", "#000075")[seq_len(K)]
names(pal) <- levels(coords$panel)

ncol_f <- if (K <= 6L) 3L else 4L
nrow_f <- ceiling(K / ncol_f)

p <- ggplot(coords, aes(x = col, y = row)) +
  ## no `panel` column -> this layer is drawn in every facet
  geom_point(data = coords[, c("col", "row")], colour = "grey85",
             size = 0.12, stroke = 0) +
  geom_point(aes(colour = panel), size = 0.18, stroke = 0) +
  scale_colour_manual(values = pal, guide = "none") +
  facet_wrap(~ panel, ncol = ncol_f) +
  scale_y_reverse() +
  coord_fixed() +
  labs(title = sprintf("NovaSeq IDC_Grade3 ds90k - %s, one panel each (k = %d)",
                       if (grepl("merged", tag)) "merged CNV clusters" else "CNV k-means clusters",
                       K)) +
  theme_void(base_size = 13) +
  theme(plot.title   = element_text(face = "bold", hjust = 0.5),
        strip.text   = element_text(face = "bold", size = 11, margin = margin(4, 0, 4, 0)),
        panel.spacing = unit(6, "pt"),
        plot.background = element_rect(fill = "white", color = NA))

xr <- diff(range(coords$col)); yr <- diff(range(coords$row))
panel_w <- 3.4
ggsave(out_png, p, width = ncol_f * panel_w,
       height = nrow_f * panel_w * yr / xr + 0.9, dpi = 300, limitsize = FALSE)
message("DONE ", out_png)
