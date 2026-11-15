## VisiumHD AP5 — 8µm epithelial object (unit of analysis for Strategy 1)
##
## Loads the spaceranger 8µm binned output, attaches the pathologist histology
## annotation (one label per 8µm bin), and keeps only the epithelial bins —
## benign + tumour, which is exactly the set InferCNV needs (benign = reference,
## tumour = observation).
##
## Output: AllEpi_8.rds.gz, read by
##   - ../infercnv/run_infercnv_8um_formal.R (as run_CNV/AP5/res8/AllEpi.rds.gz, same object)
##   - figures/high_res_st/input_data/prep/prep_visiumHD_qc.R (Fig 5 QC panel,
##     as the coarsest of the three epithelial units)
##
## Run under stcnv-r-shared (Seurat 5.1.0):
##   Rscript binning_8um_epithelial.R
suppressMessages({
  library(Seurat)
})

## paths -------------------------------------------------------------------
st_dir     <- "<COLLAB_DATA>/spatial_prostate/oxford/data/HD/H1-JQTF3M2_A1"
histo_csv  <- "<COLLAB_DATA>/spatial_prostate/oxford/data/HD/annotation/AP5_10_HD_AL.csv"
out_rds    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/binning/AllEpi_8.rds.gz"

## epithelial histology classes (everything else is Stroma*) ---------------
EPI_LEVELS <- c("Benign", "Benign (basal)",
                "Tumour Clone A", "Tumour Clone B", "Tumour - IDC-P")

## load 8µm bins + annotation ----------------------------------------------
object <- Load10X_Spatial(data.dir = st_dir, bin.size = 8)

annot8 <- read.csv(histo_csv)                       # Barcode, Histology Annotations
annot8 <- tibble::column_to_rownames(annot8, "Barcode")
object <- AddMetaData(object, annot8)

## subset to epithelium -----------------------------------------------------
epi_cells <- rownames(object@meta.data)[
  object@meta.data$Histology.Annotations %in% EPI_LEVELS]
all_epi <- subset(object, cells = epi_cells)
print(table(all_epi@meta.data$Histology.Annotations))

dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(all_epi, out_rds, compress = "gzip")
cat(sprintf("DONE: %d epithelial 8µm bins -> %s\n", ncol(all_epi), out_rds))
