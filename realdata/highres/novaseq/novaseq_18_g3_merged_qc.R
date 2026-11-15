#!/usr/bin/env Rscript
# QC of the three merged CNV clusters (novaseq_17_g3_merge_kmeans.R).
#
# Asks whether M1/M2/M3 differ for reasons other than clonal CNV:
#   depth      - nCount / nFeature per cluster (the M1-M2 split tracked depth at k=12)
#   CNV burden - mean |x - 1| per cell, from the phase-1 RDS (spot_burden, unclamped)
#   burden~depth correlation, per cluster and overall
#   composition - which res-0.3 expression clusters each M is made of
#
# Run with tuscan R (ggplot2):
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_18_g3_merged_qc.R

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

PROJ    <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir <- file.path(PROJ, "processed_data", "infercnv_IDC_Grade3_ds90k")
grp_csv <- file.path(run_dir, "IDC_Grade3_epi_groups_ds90k.csv")
lab_csv <- file.path(run_dir, "outputs", "IDC_Grade3_ds90k_cnv_merged_m3.csv")
rds     <- file.path(run_dir, "app_outputs", "novaseq_idc_grade3_ds90k.rds")
out_png <- file.path(run_dir, "outputs", "IDC_Grade3_ds90k_merged_m3_qc.png")

lb <- read.csv(lab_csv)
m  <- setNames(lb[[2]], lb[[1]])
g  <- read.csv(grp_csv, row.names = 1)[names(m), ]

message("reading phase-1 RDS for per-cell CNV burden (~4.3 GB) ...")
d <- readRDS(rds)
burden <- d$spot_burden[names(m)]
rm(d); invisible(gc())

q <- data.frame(M = factor(paste0("M", m)), nCount = g$nCount_spatial,
                nFeature = g$nFeature_spatial, burden = burden,
                expr_cluster = g$cluster_frozen)

fmt <- function(v) sprintf("%s [%s-%s]", round(median(v)),
                           round(quantile(v, .25)), round(quantile(v, .75)))
summ <- do.call(rbind, lapply(split(q, q$M), function(s) data.frame(
  n            = nrow(s),
  nCount       = fmt(s$nCount),
  nFeature     = fmt(s$nFeature),
  burden       = sprintf("%.4f", median(s$burden)),
  r_burden_dep = sprintf("%+.2f", cor(s$burden, s$nCount, method = "spearman")))))
message("\n== per merged cluster: median [IQR] ==")
print(summ)

message(sprintf("\noverall Spearman burden ~ nCount: %+.2f",
                cor(q$burden, q$nCount, method = "spearman")))
message("\n== expression-cluster composition (row %) ==")
print(round(100 * prop.table(table(q$M, q$expr_cluster), 1), 1))
message("\n== Kruskal-Wallis across M ==")
message(sprintf("  nCount p = %.3g", kruskal.test(nCount ~ M, q)$p.value))
message(sprintf("  burden p = %.3g", kruskal.test(burden ~ M, q)$p.value))

long <- rbind(data.frame(M = q$M, value = q$nCount, metric = "nCount (UMI)"),
              data.frame(M = q$M, value = q$burden, metric = "CNV burden  mean|x-1|"))
p <- ggplot(long, aes(M, value, fill = M)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c(M1 = "#E6194B", M2 = "#3CB44B", M3 = "#4363D8"),
                    guide = "none") +
  scale_y_continuous(trans = "log10") +
  labs(title = "NovaSeq IDC_Grade3 ds90k - QC of the merged CNV clusters", x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave(out_png, p, width = 9, height = 4.5, dpi = 300)
message("\nDONE ", out_png)
