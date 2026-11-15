#—— 1. 准备环境 & 加载同事的辅助函数 ——#
library(Seurat)
library(dplyr)
library(mclust)       
library(arrow)
library(semla)
source("../binning/custom_grid_binning.R")   # vendored; run this script from realdata/highres/visiumHD/rctd/
setwd('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/')

mapping_coord <- readRDS('mapping_coord_2um_8um_bin2cell_Histology.rds')
mapping_coord[1:3,1:3]

## subset epithelial ----
allowed <- c("Benign","Benign (basal)", "Tumour Clone A","Tumour - IDC-P","Tumour Clone B")

# 1) 找出所有 epithelial 的 object_id
epi_map <- mapping_coord %>%
  filter(first_type == "Epithelial")

# 2) 标记每个 object_id 是否有“不纯”注释（非 allowed 或 NA）
bad_ids <- epi_map %>%
  group_by(object_id) %>%
  summarize(has_bad = any(is.na(Histology.Annotations) | is.na(Histology.Annotations)|
                            !(Histology.Annotations %in% allowed)),
            .groups = "drop") %>%
  filter(has_bad) %>%
  pull(object_id)

# 3) 只保留既是 epithelial，且 object_id 不在 bad_ids 里的行
pure_epi_map <- epi_map %>%
  filter(!object_id %in% bad_ids,
         Histology.Annotations %in% allowed)

# 查看一下结果
pure_epi_map %>% 
  distinct(object_id) %>% 
  nrow() -> n_pure_ids
message("保留了 ", n_pure_ids, " 个纯 epithelial object_id。")

# pure_epi_map 就是你要的最纯 epithelial cells 的 mapping
write.csv(pure_epi_map,'<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/histo_decon_pure_epithelial_cell_mapping_coord.csv')

RCTD_seu <- readRDS('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/RCTD_seu.rds.gz')
pure_epi_seu <- subset(RCTD_seu,object_id%in%unique(pure_epi_map$object_id))
saveRDS(pure_epi_seu,'<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/histo_decon_pure_epi_seu.rds.gz',compress = 'gzip')
