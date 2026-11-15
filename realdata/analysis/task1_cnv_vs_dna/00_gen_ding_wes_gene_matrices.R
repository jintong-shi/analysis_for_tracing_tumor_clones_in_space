#!/usr/bin/env Rscript
# Generate per-sample WES gene-level truth matrices for Ding spatial samples that
# lack prebuilt wes_gene_* matrices (To-Do #5). Reuses Step 05's exact projection
# (project_segments_to_gene_matrices) by sourcing it with the trailing main() stripped.
#
# The Ding truth_links/wes_cnv/<sample>_wes_cnv.csv files are the FULL Ding cohort
# table (multi-sample, keyed by WES sample id). We filter to the single WES sample
# that matches each spatial sample (by patient prefix), project to gene level, and
# write wes_gene_log2ratio_matrix.csv + wes_gene_call_matrix.csv into the spatial
# sample's run_links dir, where Step 05's load_truth_matrices reads them first.

suppressPackageStartupMessages(library(data.table))

SCRIPT_DIR <- "."
DATASETS_DIR <- "<DATA_ROOT>/realdata/datasets"
GENE_POS_FILE <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
COHORT_CSV <- file.path(DATASETS_DIR,
  "Ding_nature_wes/truth_links/wes_cnv/HT268B1-Th1H3Fc2U2Z1Bs1_wes_cnv.csv")

# spatial sample -> matching WES sample id (by patient prefix)
MAP <- list(
  "HT231P1-S1H3Fc2U1Z1Bs1-H3Bs2-Test" = "HT231P1-S1H3Y1D1_1",
  "HT260C1-Th1K1Fc2U1Z1Bs1"            = "HT260C1-Th1K1A3D1_1-13",
  "HT263B1-S1H1Fc2U1Z1Bs1"             = "HT263B1-S1H1A3D1_1",
  "HT268B1-Th1H3Fc2U12Z1Bs1"           = "HT268B1-Th1H3Y1D1_1",
  "HT268B1-Th1H3Fc2U22Z1Bs1"           = "HT268B1-Th1H3Y1D1_1",
  "HT268B1-Th1H3Fc2U2Z1Bs1" = "HT268B1-Th1H3Y1D1_1",
  "HT271B1-S1H3Fc2U1Z1Bs1"  = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U2Z1Bs1"  = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U3Z1Bs1"  = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U4Z1Bs1"  = "HT271B1-S1H3Y1D1_1"
)

# --- source Step 05 helpers without running its main() ---
src <- readLines(file.path(SCRIPT_DIR, "04_evaluate_pseudobulk_vs_wes_metrics.R"))
src <- src[!grepl("^main\\(\\)\\s*$", src)]
tmp <- tempfile(fileext = ".R"); writeLines(src, tmp); source(tmp)

gene_ref <- build_gene_reference(GENE_POS_FILE)
cohort <- fread(COHORT_CSV)
stopifnot("sample" %in% names(cohort))

write_gene_csv <- function(mat, path) {
  df <- data.frame(gene = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
  names(df)[-1] <- colnames(mat)
  fwrite(df, path)
}

for (spatial in names(MAP)) {
  wes_id <- MAP[[spatial]]
  sub <- cohort[sample == wes_id]
  if (nrow(sub) == 0L) { cat(sprintf("[SKIP] %s: no rows for WES id %s\n", spatial, wes_id)); next }

  seg_tmp <- tempfile(fileext = ".csv"); fwrite(sub, seg_tmp)
  proj <- project_segments_to_gene_matrices(seg_tmp, gene_ref)
  if (is.null(proj)) { cat(sprintf("[FAIL] %s: projection returned NULL\n", spatial)); next }

  out_dir <- file.path(DATASETS_DIR, "Ding_nature_wes/run_links", spatial)
  write_gene_csv(proj$log2, file.path(out_dir, "wes_gene_log2ratio_matrix.csv"))
  write_gene_csv(proj$call, file.path(out_dir, "wes_gene_call_matrix.csv"))
  cat(sprintf("[OK] %s (WES %s): %d segs -> %d x %d genes; non-NA log2=%d\n",
              spatial, wes_id, nrow(sub), nrow(proj$log2), ncol(proj$log2),
              sum(!is.na(proj$log2[, 1]))))
}
cat("[Done] Ding wes_gene_* generation finished.\n")
