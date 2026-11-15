#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Extracts per-unit QC (raw nCount / nFeature) for the new supp-visiumHD panel f
# (per-unit QC violins, ported from fig5_merge/panels/p2_qc.R) from the three
# epithelial-unit Seurat objects. No aggregation here: only the per-unit values
# already stored in @meta.data. The medians / violins are computed later in
# scripts/supp_visiumHD.R off this CSV.
#
# Run under stcnv-r-shared (Seurat 5.1.0):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_visiumHD_qc.R
# =====================================================================
suppressWarnings(suppressMessages(library(Seurat)))
R  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data"
OD <- "../visiumHD"
dir.create(OD, showWarnings = FALSE, recursive = TRUE)

# epithelial units, coarsening depth 8µm bin -> single cell -> meta-cell
QC_OBJS <- list(
  "8µm bin"     = file.path(R, "binning/AllEpi_8.rds.gz"),
  "Single cell" = file.path(R, "deconvolution/processed_data/histo_decon_pure_epi_seu.rds.gz"),
  "Meta-cell"   = file.path(R, "pseudo_spot/epi_pseudo_r40_Smax20_seu.rds.gz"))

extract_qc <- function(label, path) {
  message("[f] ", label, " <- ", basename(path))
  obj <- readRDS(path)
  md  <- obj@meta.data
  nc  <- grep("^nCount",   colnames(md), value = TRUE)[1]
  nf  <- grep("^nFeature", colnames(md), value = TRUE)[1]
  out <- data.frame(rep = label, nCount = md[[nc]], nFeature = md[[nf]])
  rm(obj, md); invisible(gc()); out
}
qc <- do.call(rbind, Map(extract_qc, names(QC_OBJS), QC_OBJS))
data.table::fwrite(qc, file.path(OD, "qc_per_unit.csv"))
message(sprintf("wrote qc_per_unit.csv: %d rows across %d units", nrow(qc), length(QC_OBJS)))
print(table(qc$rep))
