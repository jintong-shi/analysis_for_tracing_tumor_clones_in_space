library(Seurat)
library(Matrix)
setwd('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/segmentation/bin2cell')

## load data from bin2cell ----
# counts <- ReadMtx(
#   mtx    = "matrix.mtx",
#   features = "features.tsv",feature.column = 1,
#   cells  = "barcodes.tsv"
# )
# 
# # 2. 建立基础 Seurat 对象，载入 meta.data
# meta <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)
# seu <- CreateSeuratObject(
#   counts = counts,
#   meta.data = meta
# )
# 
# # 3. 读空间坐标
coords <- read.csv("spatial.csv", row.names = 1)
# 
# coords_matrix <- as.matrix(coords)
# colnames(coords_matrix) <- c("spatial_1", "spatial_2")
# spatial_reduction <- CreateDimReducObject(embeddings = coords_matrix, 
#                                           key = "spatial_", 
#                                           assay = DefaultAssay(seu))
# 
# seu[["spatial"]] <- spatial_reduction
# DimPlot(seu, reduction = "spatial",group.by  = 'labels_joint_source')
# VlnPlot(seu,feature = c('nFeature_RNA','nCount_RNA','bin_count'),pt.size = 0)
# saveRDS(seu,'bin2cell_seu_counts.rds.gz',compress = 'gzip')
seu <- readRDS('bin2cell_seu_counts.rds.gz')
## RCTD ----
#remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)

seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu, verbose = T)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
seu <- RunUMAP(seu, reduction = "pca", return.model = T, dims = 1:30, verbose = T)

## ref
ref <- readRDS("<DATA_ROOT>/realdata_tumor_prediction/visium/data/sc_ref/PCa_atlas_annotated_Oct23_update.rds")
ref <- subset(ref,celltype_major_v2%in%c('Unassigned','Cycling'),invert=TRUE)

Idents(ref) <- "celltype_major_v2"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype_major_v2)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- seu[["RNA"]]$counts
seu_cells_hd <- colnames(seu[["RNA"]])

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 30)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD,'<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/RCTD_res_withoutCyclingUnassigned.rds')

