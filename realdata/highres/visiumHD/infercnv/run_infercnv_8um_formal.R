## VisiumHD AP5 — 8µm InferCNV formal run (benign-as-reference)
## Input: all epithelial bins (benign + tumour), Histology.Annotations as groups.
## Reference = benign epithelial; observation = tumour epithelial.
## Params: ref_group_names=c("Benign","Benign (basal)"), cluster_by_groups=FALSE,
##   analysis_mode="samples" (no subclustering — previous "subclusters" default
##   crashed in .single_tumor_leiden_subclustering / as.hclust.phylo, 2026-06-03),
##   denoise=TRUE, HMM=FALSE, no_plot/no_prelim_plot=TRUE, cutoff=0.1.
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(infercnv)
})

## paths -------------------------------------------------------------------
input_rds  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/run_CNV/AP5/res8/AllEpi.rds.gz"
gene_order <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
out_dir    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/run_CNV/AP5/res8/formal_run"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

## inputs ------------------------------------------------------------------
seu_obj.epi <- readRDS(input_rds)
matrix <- Seurat::GetAssayData(seu_obj.epi, layer = "counts", assay = "Spatial.008um") %>% as.matrix()
histology_anno <- seu_obj.epi@meta.data[, "Histology.Annotations"] %>% as.data.frame()
rownames(histology_anno) <- rownames(seu_obj.epi@meta.data)
colnames(histology_anno) <- "histology"

## infercnv object ---------------------------------------------------------
infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix = matrix,
                                         gene_order_file = gene_order,
                                         annotations_file = histology_anno,
                                         delim = "\t",
                                         ref_group_names = c("Benign", "Benign (basal)"),
                                         chr_exclude = c("chrM"))

infCNV <- infercnv::run(infCNV,
                        cutoff = 0.1,
                        out_dir = ".",
                        cluster_by_groups = FALSE,
                        analysis_mode = "samples",
                        denoise = TRUE,
                        HMM = FALSE,
                        no_plot = TRUE,
                        no_prelim_plot = TRUE)

file.copy("run.final.infercnv_obj", "visiumHD_AP5_8um_run.final.infercnv_obj", overwrite = TRUE)

cat("DONE: infercnv 8um formal run\n")
