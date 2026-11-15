#—— 1. 准备环境 & 加载同事的辅助函数 ——#
library(Seurat)
library(dplyr)
library(mclust)       
library(arrow)
library(semla)
source("../binning/custom_grid_binning.R")   # vendored; run this script from realdata/highres/visiumHD/rctd/
setwd('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/deconvolution/processed_data/without_cycling_unassigned/')

#—— 2. 读取 2µm 数据 & 8µm 注释，扩散到 2µm barcode 上 ——####
# 2.1 载入 2µm Seurat 对象（或直接用 Load10X_Spatial）
main.data.dir <- '<COLLAB_DATA>/spatial_prostate/oxford/data/HD/'
sample.name <- 'H1-JQTF3M2_A1'
st.dir <- file.path(paste0(main.data.dir,sample.name))

object <- Load10X_Spatial(data.dir = st.dir, bin.size = c(2, 8))
annot8 <- read.csv('<COLLAB_DATA>/spatial_prostate/oxford/data/HD/annotation/AP5_10_HD_AL.csv')
annot8 <- tibble::column_to_rownames(annot8,'Barcode')
object <- AddMetaData(object,annot8)
VlnPlot(object,features = c('nCount_Spatial.002um', 'nFeature_Spatial.002um'),pt.size = 0)
VlnPlot(object,features = c('nCount_Spatial.008um', 'nFeature_Spatial.008um'),pt.size = 0,assay = 'Spatial.008um')
SpatialDimPlot(object,)
# 2.2 读取 2/8µm spot 坐标

coord_path_2um <- file.path(st.dir,'binned_outputs','square_002um',"spatial/tissue_positions.parquet")
coord_df_2um <- LoadSpatialCoordinates(coordinatefiles = coord_path_2um,
                                   remove_spots_outside_tissue = FALSE,
                                   verbose = TRUE)
colnames(coord_df_2um)[1] <- 'barcode_2um'
nrow(coord_df_2um)

coord_path_8um <- file.path(st.dir,'binned_outputs','square_008um',"spatial/tissue_positions.parquet")
coord_df_8um <- LoadSpatialCoordinates(coordinatefiles = coord_path_8um,
                                       remove_spots_outside_tissue = FALSE,
                                       verbose = TRUE)
colnames(coord_df_8um)[1] <- 'barcode_8um'
nrow(coord_df_8um)
coord_df_8um_withHist <- merge(coord_df_8um,annot8,by.x = 'barcode_8um',by.y='Barcode',all=TRUE)
nrow(coord_df_8um_withHist)

# 2.3 mapping 2-8um coord 
bin_scale <- 4L

coord2 <- coord_df_2um %>%
  mutate(
    x8 = floor(x / bin_scale),
    y8 = floor(y / bin_scale)
  )

coord8 <- coord_df_8um_withHist %>%
  rename(x8 = x, y8 = y) %>%
  select(barcode_8um, x8, y8, Histology.Annotations)

mapping_coord <- coord2 %>%
  left_join(coord8, by = c("x8", "y8"))
mapping_coord[1:8, ]
nrow(mapping_coord)
mapping_coord <- mapping_coord %>% select(c('barcode_2um','barcode_8um','x','y','x8','y8','Histology.Annotations'))
saveRDS(mapping_coord,'mapping_coord_2um_8um_Histology.rds')

#—— 3. 读取 bin2cell + deconvolution 结果，并映射到 2µm barcode ——####
# 3.1 载入 bin2cell-meta和RCTD Seurat
mapping_coord <- readRDS('../mapping_coord_2um_8um_Histology.rds')
RCTD.seu <- readRDS('RCTD_seu.rds.gz')
RCTD.meta <- RCTD.seu@meta.data
RCTD.res <- readRDS('RCTD_res.rds')
RCTD.res.df <- RCTD.res@results$results_df
results <- RCTD.res@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD.res@cell_type_info$info[[2]] #list of cell type names

bin2cell_coord <- read.csv('<DATA_ROOT>/realdata_tumor_prediction/visiumHD/processed_data/segmentation/bin2cell/metadata_labels_joint.csv')
head(bin2cell_coord)
bin2cell_coord <- bin2cell_coord %>% select(X, labels_joint)

# 3.2 提取每个 spot 的 deconv 主要细胞类型
deconv_df <- RCTD.meta %>%
  select(object_id, first_type)%>%merge(bin2cell_coord,by.x = 'object_id',by.y = 'labels_joint',all=TRUE)
colnames(deconv_df)[3] <- 'barcode_2um' 

# 3.4 join 并写回 2µm Seurat metadata
mapping_coord <- mapping_coord %>%
  left_join(deconv_df, by = "barcode_2um")
saveRDS(mapping_coord,'mapping_coord_2um_8um_bin2cell_Histology.rds')

norm_weights.df <- norm_weights %>% as.data.frame() %>% tibble::rownames_to_column('object_id')

mapping_coord$object_id <- as.character(mapping_coord$object_id)
mapping_coord_deconscore <- mapping_coord %>% left_join(norm_weights.df,by = 'object_id')

#—— 4. 在同一 Seurat 对象中比较 histology vs. deconvolution ——####
mapping_coord_major <- mapping_coord_deconscore %>%
  mutate(
    ## 1) deconvolution 的 first_type → 大类
    deconv_major = case_when(
      first_type %in% c("Epithelial")                          ~ "Epi",
      first_type %in% c("B-cells","T-cells","Myeloid",
                        "Plasmablasts")  ~ "Immune",
      first_type %in% c("CAFs","Endothelial","SMCs")           ~ "Stroma",
      TRUE                                                     ~ "Other"
    ),
    ## 2) 组织学注释 Histology.Annotations → 大类
    histo_major = case_when(
      grepl("Lymphocytes|TLS", Histology.Annotations, ignore.case = TRUE) ~ "Immune",
      grepl("^Stroma", Histology.Annotations, ignore.case = TRUE)   ~ "Stroma",
      grepl("Tumour|Benign", Histology.Annotations, ignore.case = TRUE) ~ "Epi",
      TRUE                                                         ~ "Other"
    )
  )

## 5. score compare:vlnplot 不太好-----
library(tidyverse)

# 假设你的数据框叫 mapping_coord_major
df <- mapping_coord_major

# 1. 计算三类大类得分
df2 <- df %>%
  mutate(
    immune_prop  = rowSums(select(., `B-cells`, `T-cells`, Myeloid, Plasmablasts), na.rm = TRUE),
    stroma_prop  = rowSums(select(., CAFs, Endothelial, SMCs),        na.rm = TRUE),
    epi_prop     = Epithelial
  ) 

library(tidyverse)
# 1. pivot longer
long_df <- df2 %>%
  pivot_longer(
    cols      = c(immune_prop, stroma_prop, epi_prop),
    names_to  = "cell_group",
    values_to = "prop"
  )

# 2. Violin plot：分别为 immune/stroma/epi 三个 facet
p1 <- ggplot(long_df, aes(x = Histology.Annotations, y = prop, fill = cell_group)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ cell_group, scales = "free_y", ncol = 1) +
  labs(
    title = "Violin plots of immune, stroma, and epi proportions by histology",
    x     = "Histology",
    y     = "Proportion"
  ) +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(angle = 45, hjust = 1),
    strip.text         = element_text(face = "bold")
  )
png('cell_type_compare_annotation.png',400,400)
print(p1)
dev.off()

# 3. Boxplot：三组并排对比
p2 <- ggplot(long_df, aes(x = histo_major, y = prop, fill = cell_group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Boxplots of proportions by histology and cell group",
    x     = "Histology",
    y     = "Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.title    = element_blank()
  )

# 任选其一单独显示，或者并排/上下合并看：
# library(patchwork)
# p1 + p2 + plot_layout(ncol = 2)

# 直接打印：
print(p1)
print(p2)

# 2. Violin plot：immune+stroma vs epi in每个 histo-major
p1 <- df2 %>%
  ggplot(aes(x = Histology.Annotations, y = immstr_vs_epi, fill = histo_major)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ Histology.Annotations, scales = "free_y", ncol = 3) +
  labs(title = "各 Histology.Major 中 immune+stroma vs epi 得分分布",
       x = NULL, y = "Score") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# 3. Violin plot：immune vs stroma vs epi 在每个 histo-major
#    用 position_dodge 展示三组并排对比
long3 <- df2 %>%
  select(Histology.Annotations, immune_prop, stroma_prop, epi_prop) %>%
  pivot_longer(
    cols = c(immune_prop, stroma_prop, epi_prop),
    names_to  = "group3",
    values_to = "score3"
  )

p2 <- ggplot(long3, aes(x = Histology.Annotations, y = score3, fill = group3)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
  labs(title = "各 Histology.Major 中 immune vs stroma vs epi 占比分布",
       x = NULL, y = "Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

# 如果想一页上看两张图：
library(patchwork)
p1 + p2 + plot_layout(ncol = 1)


## 6. test epi accuracy ----
library(dplyr)
library(caret) 
# 假设你已经有了 mapping_coord2，包含 deconv_major 和 histo_major 两列
# 先构造二值标签
mapping_coord_major.clean <- mapping_coord_major %>% 
  drop_na(Histology.Annotations, object_id, first_type)
mapping_coord_major.clean <- mapping_coord_major.clean %>% dplyr::filter(first_type!='Unknown')

mapping_epi <- mapping_coord_major.clean %>% 
  mutate(
    pred_epi  = (deconv_major == "Epi"),       # 预测是否上皮
    truth_epi = (histo_major  == "Epi")        # 真实是否上皮
  )

#---- 简单的混淆矩阵 ----
cm <- table(
  Truth  = mapping_epi$truth_epi,
  Predict= mapping_epi$pred_epi
)
print(cm)

#---- 用 caret 计算常用指标 ----
# caret 要求 factor 且把正例放在 “positive” 参数里
mapping_epi <- mapping_epi %>%
  mutate(
    truth_epi  = factor(truth_epi,  levels=c(FALSE, TRUE)),
    pred_epi   = factor(pred_epi,   levels=c(FALSE, TRUE))
  )
conf <- confusionMatrix(
  data    = mapping_epi$pred_epi,
  reference = mapping_epi$truth_epi,
  positive = "TRUE"
)
print(conf)

#—— 7. 可视化各celltype在各位置占比 ——#####
library(ggplot2)
#mapping_coord.filter <- mapping_coord %>% dplyr::filter(first_type %in% c('CAFs','Endothelial','B-cells','SMCs',''))
conf_mat <- table(mapping_coord$Histology.Annotations, mapping_coord$first_type)

df_plot <- as.data.frame(conf_mat) %>%
  mutate(prop = Freq / tapply(Freq, Var1, sum)[Var1])

p <- ggplot(df_plot, aes(x = Var2, y = Var1, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Deconvolution label",
    y = "Histology label (8µm)",
    fill = "Row‑normalized\nproportion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf('cell_type_proportion_annotation.pdf',width = 8,height = 5)
print(p)
dev.off()

## subset epithelial -----
all_epithelial_cells <- subset(RCTD.seu,first_type=='Epithelial')
saveRDS(all_epithelial_cells,'epithelial_cells_seu.rds.gz',compress = 'gzip')
all_epithelial_cells.id <- colnames(all_epithelial_cells)
write.table(all_epithelial_cells.id,'epithelial_cells_id.txt',quote = FALSE)
