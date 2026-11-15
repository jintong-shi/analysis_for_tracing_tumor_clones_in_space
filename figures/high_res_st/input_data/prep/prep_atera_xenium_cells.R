#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Extracts the minimal raw per-cell fields the supp-atera panels consume from
# the 6.0 GB Xenium Seurat object into a compact CSV. No aggregation is done
# here: only per-cell values that already live in the object (tissue coords,
# UMAP embedding, nCount/nFeature, log-normalized EPCAM, cluster, cell type).
# All plotting/derivation happens later in scripts/supp_atera.R off this CSV.
#
# Mirrors the field set of the original supp_high_res_merge cache_df.rds.
# Run under stcnv-r-shared (Seurat v5.1.0):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_atera_xenium_cells.R
# =====================================================================
suppressWarnings(suppressMessages(library(Seurat)))

OBJ <- "<DATA_ROOT>/high_res_ST/atera/processed_data/clustering/xenium_obj.rds.gz"
OUT <- "../atera/xenium_supp_cells.csv"
stopifnot(file.exists(OBJ))

message("reading 6.0 GB gz Seurat object (slow) ...")
seu   <- readRDS(OBJ)
cells <- colnames(seu)

tc  <- GetTissueCoordinates(seu); rownames(tc) <- tc$cell; tc <- tc[cells, ]
emb <- Embeddings(seu, "umap")[cells, ]
md  <- seu@meta.data[cells, ]
vd  <- FetchData(seu, vars = c("EPCAM", "major_celltype"), layer = "data")
stopifnot(identical(rownames(vd), cells),
          "nFeature_Xenium" %in% colnames(md))

df <- data.frame(
  cell     = cells,
  x        = tc$x, y = tc$y,
  umap1    = emb[, 1], umap2 = emb[, 2],
  nCount   = md$nCount_Xenium,
  nFeature = md$nFeature_Xenium,
  EPCAM    = vd$EPCAM,
  res01    = as.character(md$Xenium_snn_res.0.1),
  celltype = as.character(md$major_celltype),
  stringsAsFactors = FALSE)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(df, OUT)
message(sprintf("wrote %s  (%d cells x %d cols)", OUT, nrow(df), ncol(df)))
