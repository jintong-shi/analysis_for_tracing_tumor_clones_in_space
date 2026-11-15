#!/usr/bin/env Rscript
# One-time builder for input_data/novaseq/ (Fig 5 panels g + h).
#
# NovaSeq X spatial, IDC_Grade3, ds90k inferCNV run (90,000 cancer cells + 476
# benign reference cells). Two extracts:
#   novaseq_g3_merged_m3_cnv_input.rds  copied verbatim from the novaseq arm; it
#     already holds mat_ref (476 x genes, unbinned), mat_bin (obs rows averaged in
#     bins of 12 within each merged CNV cluster), bin_lab, chr_genes.
#   novaseq_g3_coords_labels.csv        barcode, x, y, group — group is "Normal"
#     for the reference cells and C1/C2/C3 for the merged CNV clusters.
#
# Not a figure script; kept for provenance. Run once, from anywhere:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_novaseq_extracts.R

suppressWarnings(suppressMessages(library(data.table)))

NOVA <- "<DATA_ROOT>/high_res_ST/novaseq"
RUN  <- file.path(NOVA, "processed_data", "infercnv_IDC_Grade3_ds90k")

args <- commandArgs(trailingOnly = FALSE)
SELF <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
OUT  <- file.path(dirname(dirname(SELF)), "novaseq")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. CNV matrix bundle (already binned by the novaseq prep) ----------------
src <- file.path(RUN, "outputs", "IDC_Grade3_ds90k_allcells_merged_m3_heatmap_input.rds")
dst <- file.path(OUT, "novaseq_g3_merged_m3_cnv_input.rds")
stopifnot(file.exists(src))
file.copy(src, dst, overwrite = TRUE)
d <- readRDS(dst)
message(sprintf("cnv input: ref=%d obs=%d binned rows=%d genes=%d labels=%s",
                d$n_ref, d$n_obs, nrow(d$mat_bin), ncol(d$mat_bin),
                paste(d$lab_levels, collapse = ",")))

## ---- 2. coordinates + group label --------------------------------------------
g  <- fread(file.path(RUN, "IDC_Grade3_epi_groups_ds90k.csv"))
setnames(g, 1, "barcode")
m  <- fread(file.path(RUN, "outputs", "IDC_Grade3_ds90k_cnv_merged_m3.csv"))
lab <- setNames(paste0("C", m$merged), m$barcode)

co <- g[, .(barcode, x, y, group)]
co[, group := ifelse(group == "normal", "Normal", lab[barcode])]
if (anyNA(co$group))
  stop(sprintf("%d cells have neither a merged CNV cluster nor the normal label",
               sum(is.na(co$group))), call. = FALSE)
fwrite(co, file.path(OUT, "novaseq_g3_coords_labels.csv"))
message("coords/labels:"); print(table(co$group))
message("DONE -> ", OUT)
