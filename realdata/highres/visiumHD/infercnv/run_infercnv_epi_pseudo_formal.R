## VisiumHD AP5 — Strategy 3 (epi pseudo-spot, region-growing r40_Smax20)
## InferCNV formal run (benign-as-reference, tumour-as-observation)
##
## Unit of analysis = epithelial meta-cells (pseudo-spots) built by region
## growing (../pseudo_spot/region_growing.R -> generate_pseudo_spot_seu.R). Purity
## is precomputed in the object:
##   meta.data$status    "Pure" / "Mixed"
##   meta.data$pure_type "Benign" / "Tumor"  (NA for Mixed spots)
## Mixed spots are dropped (filter status == "Pure"):
##   - pure_type "Benign" -> reference
##   - pure_type "Tumor"  -> observation
##
## Params mirror res8/formal_run and the Strategy-2 formal run:
##   ref_group_names="Benign", cluster_by_groups=FALSE, analysis_mode="samples"
##   (no subclustering), denoise=TRUE, HMM=FALSE, cutoff=0.1, chr_exclude="chrM".
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(infercnv)
})

## paths -------------------------------------------------------------------
seu_rds    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/pseudo_spot/epi_pseudo_r40_Smax20_seu.rds.gz"
gene_order <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
out_dir    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/run_CNV/AP5/epi_pseudo_spot/r40_Smax20/formal_run"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

## subset to pure spots, build annotation ----------------------------------
seu <- readRDS(seu_rds)
n_mixed <- sum(seu@meta.data$status == "Mixed")
pure_cells <- rownames(seu@meta.data)[seu@meta.data$status == "Pure"]
seu <- subset(seu, cells = pure_cells)
cat(sprintf("pure spots: %d (dropped mixed=%d)\n", ncol(seu), n_mixed))

anno <- data.frame(histology = seu@meta.data$pure_type,
                   row.names = colnames(seu))
cat("annotation table:\n"); print(table(anno$histology))

matrix <- Seurat::GetAssayData(seu, layer = "counts", assay = "RNA") %>% as.matrix()

## infercnv object ---------------------------------------------------------
infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix = matrix,
                                         gene_order_file  = gene_order,
                                         annotations_file = anno,
                                         ref_group_names  = c("Benign"),
                                         chr_exclude      = c("chrM"))

infCNV <- infercnv::run(infCNV,
                        cutoff            = 0.1,
                        out_dir           = ".",
                        cluster_by_groups = FALSE,
                        analysis_mode     = "samples",
                        denoise           = TRUE,
                        HMM               = FALSE)

file.copy("run.final.infercnv_obj", "visiumHD_AP5_pseudospot_r40_Smax20_run.final.infercnv_obj", overwrite = TRUE)

cat("DONE: infercnv epi pseudo-spot r40_Smax20 formal run\n")