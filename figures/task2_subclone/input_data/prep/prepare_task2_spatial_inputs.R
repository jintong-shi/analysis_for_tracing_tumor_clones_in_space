#!/usr/bin/env Rscript
# =====================================================================
# FREEZE STEP (run once): build the self-contained spatial inputs for the
# Task-2 spatial panels — Fig 3 f/g and figS_task2B f/g. Pure file copy: every
# file is taken byte-identical from the pipeline outputs, nothing is recomputed.
#
# Per sample, into ../real/spatial/<dataset>/<sample>/ :
#   positions.csv          spaceranger tissue_positions.csv (full, in_tissue kept)
#   he_lowres.png          spaceranger tissue_lowres_image.png
#   keep_barcodes.txt      HistoTumor_Benign-ref run; tumour mask = keep \ ref,
#   ref_barcodes.txt         the set difference is done in the figure scripts
#   subclone__<Tool>.csv   <run>/<tool>/subclone_spots_x_subclones.csv, renamed
#
# Best-run selection (fixed, matches the assortativity best-case table):
#   reference-based tools (InferCNV/SCEVAN/Numbat) -> HistoTumor_Benign-ref
#   reference-free tools (STARCH/CalicoST/TUSCAN)  -> AllSpots_no-ref
# A tool is copied when its subclone table exists; the figure scripts render a
# "(missing)" placeholder for the ones that do not.
#
# Verified 2026-08-04: the output is byte-identical to the folder that was
# copied by hand on 2026-07-17 (diff -rq, all four samples).
# =====================================================================

# Run this from the prep/ directory; outputs go to ../real/spatial.
DATS <- "<DATA_ROOT>/realdata_unified_by_dataset/datasets"
OUT  <- "../real/spatial"

REF_BASED <- c(InferCNV = "infercnv", SCEVAN = "SCEVAN", Numbat = "NUMBAT")
REF_FREE  <- c(STARCH = "STARCH", CalicoST = "CalicoST", TUSCAN = "TUSCAN")

# run_dir / spatial_dir differ per dataset, so both are spelled out per sample.
jobs <- list(
  list(ds = "prostate_eva_ffpe", s = "P10_P3",
       run_dir     = file.path(DATS, "prostate_eva_ffpe/run_CNV/P10_P3"),
       spatial_dir = file.path(DATS, "prostate_eva_ffpe/source_links/spaceranger_P10_P3/outs/spatial")),
  list(ds = "nature_prostate_ff", s = "H1_4",
       run_dir     = file.path(DATS, "nature_prostate_ff/run_links/H1_4"),
       spatial_dir = file.path(DATS, "nature_prostate_ff/truth_links/per_sample/H1_4/spatial")),
  list(ds = "metastatic_prostate_ffpe", s = "pt15_P2",
       run_dir     = file.path(DATS, "metastatic_prostate_ffpe/run_links/pt15_P2"),
       spatial_dir = file.path(DATS, "metastatic_prostate_ffpe/source_links/selected_oxford_orig_data/pt15_P2/spatial")),
  list(ds = "metastatic_prostate_ffpe", s = "pt10_A3",
       run_dir     = file.path(DATS, "metastatic_prostate_ffpe/run_links/pt10_A3"),
       spatial_dir = file.path(DATS, "metastatic_prostate_ffpe/source_links/selected_oxford_orig_data/pt10_A3/spatial")))

copy_to <- function(src, dst) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, dst, overwrite = TRUE)
}

for (j in jobs) {
  od <- file.path(OUT, j$ds, j$s)
  dir.create(od, showWarnings = FALSE, recursive = TRUE)

  copy_to(file.path(j$spatial_dir, "tissue_positions.csv"),   file.path(od, "positions.csv"))
  copy_to(file.path(j$spatial_dir, "tissue_lowres_image.png"), file.path(od, "he_lowres.png"))
  ref_run <- file.path(j$run_dir, "HistoTumor_Benign-ref")
  copy_to(file.path(ref_run, "keep_barcodes.txt"), file.path(od, "keep_barcodes.txt"))
  copy_to(file.path(ref_run, "ref_barcodes.txt"),  file.path(od, "ref_barcodes.txt"))

  copied <- character(0)
  for (grp in list(list(run = "HistoTumor_Benign-ref", tools = REF_BASED),
                   list(run = "AllSpots_no-ref",       tools = REF_FREE))) {
    for (lab in names(grp$tools)) {
      src <- file.path(j$run_dir, grp$run, grp$tools[[lab]],
                       "subclone_spots_x_subclones.csv")
      if (copy_to(src, file.path(od, paste0("subclone__", lab, ".csv"))))
        copied <- c(copied, lab)
    }
  }
  cat(sprintf("[freeze] %s/%s  tools: %s\n", j$ds, j$s, paste(copied, collapse = ", ")))
}
