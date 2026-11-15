## VisiumHD AP5 — epithelial meta-cells by region growing (step 2/2)
##
## Aggregates the cells of each meta-cell from step 1 (region_growing.R) into one
## pseudo-spot: summed counts, centroid coordinates as the "spatial" reduction,
## and the histology purity call InferCNV keys off.
##
## Purity: every cell carries an 8µm histology label; collapse it to
## Benign vs Tumor, then a meta-cell is "Pure" if all its annotated cells agree
## ("pure_type" = that label) and "Mixed" otherwise (pure_type = NA). The
## InferCNV runner (../infercnv/run_infercnv_epi_pseudo_formal.R) drops Mixed and uses
## Benign as reference / Tumor as observation.
##
## Output: epi_pseudo_r<r>_Smax<S>_seu.rds.gz, also read by the Fig-5 QC prep
## (figures/high_res_st/input_data/prep/prep_visiumHD_qc.R) as the finest of the
## three epithelial units.
##
## Run under stcnv-r-shared:
##   Rscript generate_pseudo_spot_seu.R
suppressMessages({
  library(Seurat)
  library(dplyr)
})

## paths / parameter --------------------------------------------------------
STEM    <- "r40_Smax20"
outdir  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/pseudo_spot"
map_csv <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/deconvolution/processed_data/histo_decon_pure_epithelial_cell_mapping_coord.csv"

## per-cell histology, Benign (basal) collapsed into Benign ------------------
epi_map <- read.csv(map_csv, row.names = 1) %>%
  filter(first_type == "Epithelial", !is.na(Histology.Annotations)) %>%
  mutate(Histology.Annotations = ifelse(
    Histology.Annotations %in% c("Benign", "Benign (basal)"),
    "Benign", Histology.Annotations)) %>%
  select(object_id, Histology.Annotations) %>%
  mutate(object_id = as.character(object_id))   # colnames() is character; dplyr won't join across types

## meta-cell membership from step 1 -----------------------------------------
epi_pseudo_cells <- readRDS(file.path(outdir, sprintf("epi_pseudo_%s.rds", STEM)))
epi_pseudo_cells <- subset(epi_pseudo_cells,
                           cells = colnames(epi_pseudo_cells)[!is.na(epi_pseudo_cells$pseudo_id)])

## aggregate counts ---------------------------------------------------------
epi_pseudo_spots <- AggregateExpression(epi_pseudo_cells, group.by = "pseudo_id",
                                        assays = "RNA", slot = "counts",
                                        return.seurat = TRUE)
mat <- epi_pseudo_spots@assays$RNA@layers$counts
epi_pseudo_spots$nCount_RNA   <- Matrix::colSums(mat)
epi_pseudo_spots$nFeature_RNA <- Matrix::colSums(mat > 0)

## centroid coordinates -> "spatial" reduction ------------------------------
coords_df <- as.data.frame(epi_pseudo_cells@reductions$spatial@cell.embeddings)
colnames(coords_df)[1:2] <- c("x", "y")
coords_df$pseudo_id <- epi_pseudo_cells$pseudo_id

centroids <- coords_df %>%
  group_by(pseudo_id) %>%
  summarize(x_centroid = mean(x), y_centroid = mean(y), .groups = "drop") %>%
  as.data.frame()
rownames(centroids) <- paste0("g", centroids$pseudo_id)

coords_mat <- as.matrix(centroids[colnames(epi_pseudo_spots),
                                  c("x_centroid", "y_centroid")])
colnames(coords_mat) <- c("spatial_1", "spatial_2")
epi_pseudo_spots[["spatial"]] <- CreateDimReducObject(
  embeddings = coords_mat, key = "spatial_",
  assay = DefaultAssay(epi_pseudo_spots))

## purity call --------------------------------------------------------------
## Joined on pseudo_id, not positionally: meta-cells whose cells are all
## unannotated do not appear in `purity` and stay NA.
purity <- data.frame(object_id = colnames(epi_pseudo_cells),
                     pseudo_id = epi_pseudo_cells$pseudo_id) %>%
  left_join(epi_map, by = "object_id") %>%
  filter(!is.na(Histology.Annotations)) %>%
  mutate(Ann2 = ifelse(Histology.Annotations == "Benign", "Benign", "Tumor")) %>%
  group_by(pseudo_id) %>%
  summarize(n_anno = n_distinct(Ann2), types = paste(unique(Ann2), collapse = "; "),
            .groups = "drop") %>%
  mutate(status    = ifelse(n_anno > 1, "Mixed", "Pure"),
         pure_type = ifelse(status == "Pure", types, NA_character_)) %>%
  as.data.frame()                               # tibbles do not carry row names
rownames(purity) <- paste0("g", purity$pseudo_id)

epi_pseudo_spots$status    <- purity[colnames(epi_pseudo_spots), "status"]
epi_pseudo_spots$pure_type <- purity[colnames(epi_pseudo_spots), "pure_type"]
print(table(epi_pseudo_spots$status, useNA = "ifany"))
print(table(epi_pseudo_spots$pure_type, useNA = "ifany"))

out_rds <- file.path(outdir, sprintf("epi_pseudo_%s_seu.rds.gz", STEM))
saveRDS(epi_pseudo_spots, out_rds, compress = "gzip")
cat(sprintf("DONE: %d pseudo-spots -> %s\n", ncol(epi_pseudo_spots), out_rds))