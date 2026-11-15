## VisiumHD AP5/pt10 — regular-Visium ground-truth InferCNV run
## (benign-as-reference, tumour-as-observation)
##
## The GT row of the Fig 5 b/c strategy comparison, and the benign-spot
## background of panel d: the same pt10_AP5 tissue profiled with regular
## 55 µm Visium (metastatic_prostate_ffpe dataset), run through inferCNV with
## the pathologist three-category histology (Tumor / Normal / Others).
## Reference = "Normal" spots, observation = "Tumor" spots; "Others" is dropped
## by the subset.
##
## Params are NOT the same as the three HD-strategy scripts — this is the run
## as executed: analysis_mode="cells" (spot-level, no subclustering),
## num_threads=4, no chr_exclude, plots kept (writes infercnv.png). Do not
## harmonise them with run_infercnv_8um_formal.R or the GT row stops matching
## the frozen object.
##
## Output: <out_dir>/run.final.infercnv_obj, frozen as
##   figures/high_res_st/input_data/visiumHD/infercnv/visium_pt10_AP5_GT_run.final.infercnv_obj
## Originally lived at processed_data/run_CNV/visium/run_cnv.R.
library(infercnv)
library(dplyr)
library(Seurat)
library(readr)

## paths -------------------------------------------------------------------
DATS       <- "<DATA_ROOT>/realdata_unified_by_dataset/datasets/metastatic_prostate_ffpe"
visium_dir <- file.path(DATS, "source_links/selected_oxford_orig_data/pt10_AP5")
histo_csv  <- file.path(DATS, "truth_links/histology/pt10_AP5_three_categories.csv")
gene_pos_file <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
out_dir    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/run_CNV/AP5/visium/"

## inputs ------------------------------------------------------------------
AP5_seu <- Seurat::Load10X_Spatial(visium_dir)
AP5_hist <- read_csv(histo_csv)
AP5_seu <- AddMetaData(AP5_seu, AP5_hist %>% tibble::column_to_rownames('Barcode'))
AP5_seu_epi <- subset(AP5_seu, histology %in% c('Normal', 'Tumor'))

# create the infercnv object
counts <- GetAssayData(AP5_seu_epi, assay = 'Spatial', layer = "counts")
annot_df <- data.frame(row.names = rownames(AP5_seu_epi@meta.data), group = AP5_seu_epi@meta.data[, "histology"])

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = annot_df,
  gene_order_file   = gene_pos_file,
  delim             = "\t",
  ref_group_names   = 'Normal'
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = FALSE,
  num_threads = 4,
  analysis_mode = 'cells',
  HMM = FALSE,
  denoise = TRUE
)

file.copy(file.path(out_dir, "run.final.infercnv_obj"),
          file.path(out_dir, "visium_pt10_AP5_GT_run.final.infercnv_obj"), overwrite = TRUE)
