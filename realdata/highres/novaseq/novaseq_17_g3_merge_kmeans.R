#!/usr/bin/env Rscript
# Merge the k=12 CNV k-means clusters into the final clusters.
#
# Grouping is the k=3 cut of the CNV-profile dendrogram in
# novaseq_16_g3_kmeans_merge_stats.R (average linkage on 1 - correlation of the
# mean profiles), kept only because each group is also spatially contiguous:
#
#   M1 = C1 C2 C6 C7 C10 C11 C12   left / lower tissue mass, 52,984 cells
#        within-group CNV r 0.86-0.98; spatial cosine up to 0.79, centroids
#        all within 0.22 of the tissue diagonal
#   M2 = C4 C5 C8 C9               upper-right lobe and band, 32,609 cells
#        within-group CNV r 0.91-0.96; centroids within 0.16
#   M3 = C3                        4,407 cells, the CNV outlier (r 0.68-0.88 to
#        everything else), spatially diffuse -- left on its own
#
# Run with tuscan R:
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_17_g3_merge_kmeans.R

run_dir <- paste0("<DATA_ROOT>/high_res_ST/",
                  "novaseq/processed_data/infercnv_IDC_Grade3_ds90k")
in_csv  <- file.path(run_dir, "outputs", "IDC_Grade3_ds90k_cnv_kmeans_k12.csv")
out_csv <- file.path(run_dir, "outputs", "IDC_Grade3_ds90k_cnv_merged_m3.csv")

merge_map <- c("1" = 1, "2" = 1, "6" = 1, "7" = 1, "10" = 1, "11" = 1, "12" = 1,
               "4" = 2, "5" = 2, "8" = 2, "9" = 2,
               "3" = 3)

lb <- read.csv(in_csv)
lb$merged <- unname(merge_map[as.character(lb[[2]])])
stopifnot(!anyNA(lb$merged))
print(table(k12 = lb[[2]], merged = lb$merged))
print(table(lb$merged))

write.csv(data.frame(barcode = lb[[1]], merged = lb$merged), out_csv, row.names = FALSE)
message("DONE ", out_csv)
