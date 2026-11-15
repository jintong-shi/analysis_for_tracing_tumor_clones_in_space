# Source the repo's palette before setwd() moves us to the data root.
# Run this script from its own directory.
source('../../../../figures/utils/color.R')
setwd('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/')
library(Seurat)
library(spacexr)
RCTD.res <- readRDS('RCTD_res.rds')

results <- RCTD.res@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD.res@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD.res@spatialRNA

## plotting ----
resultsdir <- '<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
setwd(resultsdir)
# make the plots 
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

## add to seu ----
seu <- readRDS('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/segmentation/bin2cell/bin2cell_seu_counts.rds.gz')
seu <- AddMetaData(seu, metadata = RCTD.res@results$results_df)
seu$first_type[is.na(seu$first_type)] <- "Unknown"
table(seu$first_type)

pdf('cell_type_spatialdimplot.pdf',width = 9,height = 7)
DimPlot(seu,group.by = 'first_type',reduction = 'spatial',cols = .cluster_cols)
dev.off()

saveRDS(seu,'RCTD_seu.rds.gz',compress = 'gzip')

