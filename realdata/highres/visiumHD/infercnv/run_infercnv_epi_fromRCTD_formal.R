## VisiumHD AP5 — Strategy 2 (epithelial_fromRCTD) InferCNV formal run
## (benign-as-reference, tumour-as-observation)
##
## Unit of analysis = bin2cell-segmented epithelial cells (RCTD pure
## epithelial). The Strategy-2 Seurat object has no Histology.Annotations
## column, so each cell's histology is derived from its constituent 2um-bin
## annotations in histo_decon_pure_epithelial_cell_id.txt, aggregated per
## object_id:
##   - pure benign  (all bins Benign / Benign (basal))      -> reference
##   - pure tumour  (all bins Tumour Clone A / Clone B)      -> observation
##   - mixed        (bins span both benign AND tumour)       -> DROPPED
##   sub-label = dominant (most frequent) bin annotation within the cell.
##
## Params mirror res8/formal_run/run_infercnv_8um_formal.R:
##   ref_group_names=c("Benign","Benign (basal)"), cluster_by_groups=FALSE,
##   analysis_mode="samples" (no subclustering), denoise=TRUE, HMM=FALSE,
##   no_plot/no_prelim_plot=TRUE, cutoff=0.1, chr_exclude="chrM".
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(infercnv)
})

## paths -------------------------------------------------------------------
seu_rds    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/histo_decon_pure_epi_seu.rds.gz"
meta_txt   <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/histo_decon_pure_epithelial_cell_id.txt"
gene_order <- "<REF_ROOT>/gencode_v38_gene_pos.txt"
out_dir    <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/run_CNV/AP5/epithelial_fromRCTD/formal_run"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

## classify cells from per-bin histology -----------------------------------
benign_lbls <- c("Benign", "Benign (basal)")
tumour_lbls <- c("Tumour Clone A", "Tumour Clone B")

md <- read.csv(meta_txt, stringsAsFactors = FALSE)
md$coarse <- ifelse(md$Histology.Annotations %in% benign_lbls, "benign",
              ifelse(md$Histology.Annotations %in% tumour_lbls, "tumour", NA))
stopifnot(!any(is.na(md$coarse)))

per_cell <- md %>% group_by(object_id) %>% summarise(
  has_benign = any(coarse == "benign"),
  has_tumour = any(coarse == "tumour"),
  dom_label  = names(sort(table(Histology.Annotations), decreasing = TRUE))[1],
  .groups    = "drop")
per_cell$class <- with(per_cell, ifelse(has_benign & has_tumour, "mixed",
                                 ifelse(has_benign, "pure_benign", "pure_tumour")))
pure <- per_cell %>% filter(class != "mixed")
cat(sprintf("cells: pure_benign=%d pure_tumour=%d (dropped mixed=%d)\n",
            sum(pure$class == "pure_benign"), sum(pure$class == "pure_tumour"),
            sum(per_cell$class == "mixed")))

## subset Seurat to pure cells, build annotation ---------------------------
seu <- readRDS(seu_rds)
ids <- as.character(pure$object_id)
ids <- ids[ids %in% colnames(seu)]
seu <- subset(seu, cells = ids)

anno <- data.frame(
  histology = pure$dom_label[match(colnames(seu), as.character(pure$object_id))],
  row.names = colnames(seu))
cat("annotation table:\n"); print(table(anno$histology))

matrix <- Seurat::GetAssayData(seu, layer = "counts", assay = "RNA") %>% as.matrix()

## infercnv object ---------------------------------------------------------
infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix = matrix,
                                         gene_order_file  = gene_order,
                                         annotations_file = anno,
                                         ref_group_names  = c("Benign", "Benign (basal)"),
                                         chr_exclude      = c("chrM"))

infCNV <- infercnv::run(infCNV,
                        cutoff            = 0.1,
                        out_dir           = ".",
                        cluster_by_groups = FALSE,
                        analysis_mode     = "samples",
                        denoise           = TRUE,
                        HMM               = FALSE,
                        no_plot           = TRUE,
                        no_prelim_plot    = TRUE)

file.copy("run.final.infercnv_obj", "visiumHD_AP5_epithelial_fromRCTD_run.final.infercnv_obj", overwrite = TRUE)

cat("DONE: infercnv epithelial_fromRCTD formal run\n")
