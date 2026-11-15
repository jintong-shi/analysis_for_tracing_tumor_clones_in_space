#!/usr/bin/env Rscript
# Spatial assortativity of the merged CNV clusters (M1/M2/M3) of the ds90k run.
#
# Uses the benchmark's own implementation rather than a reimplementation: the
# reference script is sourced with its trailing main() call stripped, so
# build_knn_edges(), assortativity_nominal_undirected() and
# permute_assortativity() are byte-identical to
#   opensource/benchmark_spatialCNV/realdata/analysis/task2_subclone_vs_spatial/
#     07_evaluate_subclone_spatial_assortativity.R
# and its defaults (k = 6 neighbours, 1000 permutations, seed 123) are reused.
#
# Newman nominal assortativity on an undirected kNN graph over cell centroids;
# null = 1000 label shuffles on the fixed graph; two-sided p as in the reference.
# The reference restricts to a tumour mask; here every cell in the run is already
# a cancer cell, so the whole set is used.
#
# Run with tuscan R (has RANN + data.table):
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_19_g3_merged_assortativity.R

suppressPackageStartupMessages({ library(data.table) })

.self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
REF <- file.path(dirname(.self), "..", "..", "analysis",
                 "task2_subclone_vs_spatial", "07_evaluate_subclone_spatial_assortativity.R")
src <- readLines(REF, warn = FALSE)
eval(parse(text = src[!grepl("^\\s*main\\(\\)\\s*$", src)]), envir = globalenv())

PROJ    <- "<DATA_ROOT>/high_res_ST/novaseq"
run_dir <- file.path(PROJ, "processed_data", "infercnv_IDC_Grade3_ds90k")
grp_csv <- file.path(run_dir, "IDC_Grade3_epi_groups_ds90k.csv")
out_dir <- file.path(run_dir, "outputs")
out_csv <- file.path(out_dir, "IDC_Grade3_ds90k_spatial_assortativity.csv")

labelings <- c(merged_m3 = "IDC_Grade3_ds90k_cnv_merged_m3.csv",
               kmeans_k6 = "IDC_Grade3_ds90k_cnv_kmeans_k6.csv",
               kmeans_k12 = "IDC_Grade3_ds90k_cnv_kmeans_k12.csv")

g <- read.csv(grp_csv, row.names = 1)

rows <- lapply(names(labelings), function(nm) {
  lb  <- read.csv(file.path(out_dir, labelings[[nm]]))
  lab <- setNames(as.character(lb[[2]]), lb[[1]])
  bc  <- intersect(rownames(g), names(lab))

  dt <- data.table(spot_id = bc, x = g[bc, "x"], y = g[bc, "y"],
                   subclone = lab[bc])
  setorder(dt, spot_id)                       # same ordering as the reference

  set.seed(DEFAULT_SEED)
  edges <- build_knn_edges(dt[, .(x, y)], k = DEFAULT_K)
  labels_int <- as.integer(factor(dt$subclone, levels = sort(unique(dt$subclone))))
  o <- permute_assortativity(labels_int, edges, n_perm = DEFAULT_N_PERM)

  message(sprintf("%-11s n=%d  clusters=%d  edges=%d  r=%.4f  null=%.4f+-%.4f  z=%.1f  p=%.3g",
                  nm, nrow(dt), max(labels_int), nrow(edges),
                  o$obs, o$null_mean, o$null_sd, o$z, o$p_two_sided))

  data.table(labeling = nm, n_cells = nrow(dt), n_clusters = max(labels_int),
             n_edges = nrow(edges), k_neighbors = DEFAULT_K,
             n_permutation = DEFAULT_N_PERM,
             assortativity_obs = o$obs, null_mean = o$null_mean,
             null_sd = o$null_sd, z_score = o$z, p_value_two_sided = o$p_two_sided)
})

res <- rbindlist(rows)
fwrite(res, out_csv)
message("\nDONE ", out_csv)
