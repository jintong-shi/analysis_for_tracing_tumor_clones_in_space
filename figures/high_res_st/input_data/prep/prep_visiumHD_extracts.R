#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Extracts minimal raw per-element fields for supp-visiumHD panels c/d/e from
# the deconvolution Seurat objects / mapping RDS. No aggregation here: only the
# per-cell/per-bin values already stored in the objects. All counting / confusion
# / normalization happens later in scripts/supp_visiumHD.R off these CSVs.
#
# Run under stcnv-r-shared (Seurat 5.1.0):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_visiumHD_extracts.R
# =====================================================================
suppressWarnings(suppressMessages(library(Seurat)))
R  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data"
OD <- "../visiumHD"
dir.create(OD, showWarnings = FALSE, recursive = TRUE)

## c: RCTD single-cell deconvolution cell types (spatial embedding + first_type)
message("[c] RCTD_seu -> rctd_cells.csv")
seu <- readRDS(file.path(R, "deconvolution/RCTD_seu.rds.gz"))
emb <- seu@reductions$spatial@cell.embeddings
ft  <- as.character(seu@meta.data$first_type)
data.table::fwrite(data.frame(x = emb[, 1], y = emb[, 2], first_type = ft),
                   file.path(OD, "rctd_cells.csv"))
rm(seu, emb, ft); invisible(gc())

## d: pure-epi single-cell spatial embedding (cell -> x,y); counts done in script
message("[d] histo_decon_pure_epi_seu -> epi_cell_embedding.csv")
s   <- readRDS(file.path(R, "deconvolution/processed_data/histo_decon_pure_epi_seu.rds.gz"))
e2  <- s@reductions$spatial@cell.embeddings
data.table::fwrite(data.frame(cell = rownames(e2), x = e2[, 1], y = e2[, 2]),
                   file.path(OD, "epi_cell_embedding.csv"))
rm(s, e2); invisible(gc())

## e: deconv label vs 8um histology label per mapped bin (confusion done in script)
message("[e] mapping_coord -> deconv_vs_histology.csv")
mc <- readRDS(file.path(R, "deconvolution/mapping_coord_2um_8um_bin2cell_Histology.rds"))
data.table::fwrite(data.frame(Histology = as.character(mc$Histology.Annotations),
                              Deconv    = as.character(mc$first_type)),
                   file.path(OD, "deconv_vs_histology.csv"))
message("DONE visiumHD extracts")
