#!/usr/bin/env Rscript
# Build Standardized_Inputs for one replicate: per-tool counts, references,
# annotations and BAF matrices derived from the replicate's synthetic/ data.
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(Matrix)
  library(tibble)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript distribute_inputs.R <Replicate_Dir>")
rep_dir <- args[1]

# --- Paths ---
expr_file    <- file.path(rep_dir,'synthetic',"spatial_data", "st-expression.parquet")
meta_file    <- file.path(rep_dir,'synthetic',"spatial_data", "st-meta.parquet")
annot_file   <- file.path(rep_dir, 'synthetic',"spatial_data", "st-annotation.parquet")
sc_meta_file <- file.path(rep_dir, 'synthetic',"single_cell_data", "sc-meta_data.tsv")
map_file     <- file.path(rep_dir,'synthetic', "spatial_data", "cell_by_spot.txt")
allele_file  <- file.path(rep_dir, "numbat_allele_counts.tsv.gz")
out_root     <- file.path(rep_dir, "Standardized_Inputs")

message(paste0("\n=== Distributing Inputs: ", basename(rep_dir), " ==="))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Basic Data ---
if (!file.exists(annot_file)) stop(paste("Missing:", annot_file))
if (!file.exists(sc_meta_file)) stop(paste("Missing:", sc_meta_file))
if (!file.exists(map_file)) stop(paste("Missing:", map_file))

message("-> Loading Annotation & Single Cell Meta...")
# Spot Annotation (contains 'index' and 'annotation' -> benign/aberrant)
df_annot <- read_parquet(annot_file)
colnames(df_annot)[1] <- 'cell_id'
# Single Cell Metadata (contains '...1' as cell_id and 'cell_type')
df_sc_meta <- read_tsv(sc_meta_file, show_col_types = FALSE)
# Handle rownames column if it exists as '...1'
colnames(df_sc_meta)[1] <- 'cell_id'

# --- 2. Calculate Proportions from cell_by_spot.txt ---
message("-> Calculating Spot Proportions...")
lines <- readLines(map_file)

# Parse text file: "Spot_X: Cell_A, Cell_B"
parsed_list <- lapply(lines, function(line) {
  parts <- str_split(line, "\\s*[:]+\\s*", simplify = TRUE)
  if (length(parts) < 2) return(NULL)
  
  spot_id <- str_trim(parts[1])
  # Your snippet suggests adding "Cell_" prefix if raw IDs are numbers, 
  # but checking sc_meta first is safer. We assume sc_meta IDs match cell_by_spot IDs.
  cell_ids <- str_split(parts[2], ",")[[1]] %>% str_trim()
  
  # Check format matching based on your previous code snippet "Cell_" prefix logic
  # If sc_meta IDs start with "Cell_", ensure we match that format.
  if (any(grepl("^Cell_", df_sc_meta$cell_id[1:5])) && !any(grepl("^Cell_", cell_ids[1:5]))) {
    cell_ids <- paste0("Cell_", cell_ids)
  }
  
  return(list(spot = spot_id, cells = cell_ids))
})

# Filter nulls and convert to useful format
parsed_list <- Filter(Negate(is.null), parsed_list)

# Map Cells to Types
# Create a lookup vector for speed
sc_lookup <- setNames(df_sc_meta$cell_type, df_sc_meta$cell_id)

prop_rows <- lapply(parsed_list, function(item) {
  spot <- item$spot
  cells <- item$cells
  
  # Get types for these cells
  types <- sc_lookup[cells]
  types <- types[!is.na(types)] # Remove unmapped
  
  if (length(types) == 0) return(NULL)
  
  # Calculate proportions
  prop_table <- table(types) / length(types)
  as.data.frame(prop_table) %>%
    mutate(spot = spot) %>%
    dplyr::rename(cell_type = types, proportion = Freq)
})

# Combine and Pivot
df_props_long <- bind_rows(prop_rows)
df_props_wide <- df_props_long %>%
  pivot_wider(names_from = cell_type, values_from = proportion, values_fill = 0) %>%
  dplyr::rename(cell_id = spot)

# Rename columns to *_prop to match your convention if needed, or identify Tumor column
# We need to find the "Tumor" proportion. 
# Usually named "Tumor" or we check which column corresponds to aberrant behavior.
# We will sum up any column containing "Tumor" (case insensitive) just to be safe.
tumor_cols <- colnames(df_props_wide)[grepl("Tumor", colnames(df_props_wide), ignore.case = TRUE)]

if (length(tumor_cols) > 0) {
  df_props_wide$Tumor_prop <- rowSums(df_props_wide[, tumor_cols, drop=FALSE])
} else {
  # If no "Tumor" column, use the annotation file to guide us later, 
  # or assume anything not "Normal" is tumor? 
  # For now, initialize 0 if missing.
  warning("No 'Tumor' cell type found in SC metadata. Tumor_prop set to 0 (relying on 'aberrant' label).")
  df_props_wide$Tumor_prop <- 0
}

if ("Hepatocyte" %in% colnames(df_props_wide)) {
  df_props_wide$Hepatocyte_prop <- as.numeric(df_props_wide$Hepatocyte)
} else {
  warning("No 'Hepatocyte' cell type found in SC metadata. Hepatocyte_prop set to 0.")
  df_props_wide$Hepatocyte_prop <- 0
}

# --- 3. Merge and Define Groups ---
message("-> Merging Metadata...")
df_final <- df_annot %>%
  left_join(df_props_wide, by = "cell_id") %>%
  mutate(
    Tumor_prop = coalesce(Tumor_prop, 0),
    Hepatocyte_prop = coalesce(Hepatocyte_prop, 0)
  ) # Fill NAs with 0

# Define Groups
# Ref   = Hepatocyte>0.5 & Tumor_prop==0 
# Tumor = Tumor_prop > 0.2
df_final <- df_final %>%
  mutate(group = case_when(
    Hepatocyte_prop>0.5 & Tumor_prop==0 ~ "Normal",
    Tumor_prop > 0.2 ~ "Tumor",
    TRUE ~ "Ignore"
  ))
write.csv(df_final,paste0(out_root,'/Metadata_cellprop.csv'),row.names = FALSE)
stats <- table(df_final$group)
message(paste("-> Group Stats:", paste(names(stats), stats, sep="=", collapse=", ")))

# --- 4. Load Expression Matrix ---
message("-> Loading Expression...")
df_expr <- read_parquet(expr_file)

if ("index" %in% colnames(df_expr)) {
  gene_names <- df_expr$index
  mat_expr <- as.matrix(df_expr %>% select(-index))
  rownames(mat_expr) <- gene_names
} else {
  mat_expr <- as.matrix(df_expr)
}
mat_sparse <- as(mat_expr, "dgCMatrix")
rm(df_expr, mat_expr); gc()

# --- 5. Align and Export ---
common_cells <- intersect(colnames(mat_sparse), df_final$cell_id)
mat_sparse <- mat_sparse[, common_cells]
df_final   <- df_final %>% filter(cell_id %in% common_cells) %>% arrange(match(cell_id, colnames(mat_sparse)))

# A. Ref-mode ----
keep_mask <- df_final$group %in% c("Normal", "Tumor")
keep_mask_cell <- df_final[keep_mask,]%>%dplyr::pull('cell_id')
mat_ref   <- mat_sparse[, keep_mask]
df_ref    <- df_final[keep_mask, ]

# InferCNV (ref-mode) ----
dir_icnv <- file.path(out_root, "InferCNV")
dir.create(dir_icnv, recursive=TRUE, showWarnings=FALSE)
saveRDS(mat_ref, file.path(dir_icnv, "counts.rds"))
write.table(df_ref %>% select(cell_id, group)%>%tibble::column_to_rownames('cell_id'), 
            file.path(dir_icnv, "annotations.txt"), 
            sep="\t", quote=FALSE)

# CopyKAT ----
dir_ck <- file.path(out_root, "CopyKAT")
dir.create(dir_ck, recursive=TRUE, showWarnings=FALSE)

# 1. Save Counts (RDS format)
# We reuse mat_ref (the subset containing only Tumor + Normal cells)
saveRDS(mat_ref, file.path(dir_ck, "counts.rds"))

# 2. Save Reference Barcodes
# Extract cell IDs where group is "Normal" to use as 'norm.cell.names'
ref_barcodes <- df_ref %>% 
  filter(group == "Normal") %>% 
  pull(cell_id)

saveRDS(ref_barcodes, file.path(dir_ck, "references.rds"))

# Numbat ----
dir_nb <- file.path(out_root, "Numbat")
dir.create(dir_nb, recursive=TRUE, showWarnings=FALSE)
saveRDS(mat_ref, file.path(dir_nb, "counts.rds"))
write_tsv(df_ref %>% select(cell_id, group), 
          file.path(dir_nb, "annotations.tsv"))
if (file.exists(allele_file)) {
  target <- file.path(dir_nb, "alleles.tsv.gz")
  df_allele <- read_tsv(allele_file, show_col_types = FALSE)
  df_allele <- df_allele %>% 
    filter(cell %in% keep_mask_cell)
  write_tsv(df_allele,target)
}

# SCEVAN ----
dir_sc <- file.path(out_root, "SCEVAN")
dir.create(dir_sc, recursive=TRUE, showWarnings=FALSE)

# 1. Save Counts (RDS format)
# We reuse mat_ref (subset containing only Tumor + Normal cells)
saveRDS(mat_ref, file.path(dir_sc, "counts.rds"))

# 2. Save Reference Barcodes
# Extract cell IDs where group is "Normal"
ref_barcodes_sc <- df_ref %>% 
  filter(group == "Normal") %>% 
  pull(cell_id)

saveRDS(ref_barcodes_sc, file.path(dir_sc, "references.rds"))

# CONICSmat ----
dir_cm <- file.path(out_root, "CONICSmat")
dir.create(dir_cm, recursive=TRUE, showWarnings=FALSE)

# 1. Save Counts (RDS)
saveRDS(mat_ref, file.path(dir_cm, "counts.rds"))

# 2. Save Reference Barcodes
ref_barcodes_cm <- df_ref %>% 
  filter(group == "Normal") %>% 
  pull(cell_id)

saveRDS(ref_barcodes_cm, file.path(dir_cm, "references.rds"))

# InferCNVPy ----
dir_ipy <- file.path(out_root, "InferCNVPy")
dir.create(dir_ipy, recursive=TRUE, showWarnings=FALSE)

# Save Reference Barcodes (Text file, one per line)
ref_barcodes_ipy <- df_ref %>% 
  filter(group == "Normal") %>% 
  pull(cell_id)

write.table(ref_barcodes_ipy, file.path(dir_ipy, "ref_spots.txt"), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# XClone ----
dir_xc_root <- file.path(out_root, "XClone")
dir_xc_rdr  <- file.path(dir_xc_root, "RDR")
dir_xc_baf  <- file.path(dir_xc_root, "BAF")
dir_xc_anno <- file.path(dir_xc_root, "Annotation")

# Create directories
for(d in c(dir_xc_rdr, dir_xc_baf, dir_xc_anno)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# 1. Prepare RDR (Expression Data)
# Source: Raw synthetic data (assuming these variables are defined earlier in your script)
# You might need to adjust 'dir_sim_spatial' to match your 'synthetic/spatial_data' path
file.copy(file.path(expr_file), 
          file.path(dir_xc_rdr, "st-expression.parquet"), overwrite=TRUE)
if (file.exists(meta_file)) {
  file.copy(file.path(meta_file),
            file.path(dir_xc_rdr, "st-meta.parquet"), overwrite=TRUE)
}

# 2. Prepare BAF (Allele Data)
# Source: Standardized_Inputs/Allele_Matrix (Assuming this exists from previous steps)
# If not, point this to your cellSNP output folder
src_baf <- file.path(out_root, "Allele_Matrix")
if (!dir.exists(src_baf)) {
  src_baf <- file.path(rep_dir, "baf")
}
baf_files <- c("AD.mtx", "DP.mtx", "snps.tsv", "barcodes.tsv")

for (f in baf_files) {
  src_f <- file.path(src_baf, f)
  if (file.exists(src_f)) {
    file.copy(src_f, file.path(dir_xc_baf, f), overwrite=TRUE)
  } else {
    warning(paste("Missing BAF file:", src_f))
  }
}

# 3. Prepare Annotation
# We save a clean CSV specifically for XClone
df_xc_anno <- df_ref %>% select(cell_id, group)
write.csv(df_xc_anno, file.path(dir_xc_root, "annotation.csv"), row.names=FALSE, quote=FALSE)
write.csv(df_xc_anno, file.path(dir_xc_anno, "annotation.csv"), row.names=FALSE, quote=FALSE)

message(">>> XClone inputs prepared in: ", dir_xc_root)

# CalicoST ----
dir_ca_root    <- file.path(out_root, "CalicoST")
dir_ca_baf     <- file.path(dir_ca_root, "baf")
dir_ca_spatial <- file.path(dir_ca_root, "spatial")
for (d in c(dir_ca_root, dir_ca_baf, dir_ca_spatial)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

src_baf <- file.path(out_root, "Allele_Matrix")
if (!dir.exists(src_baf)) {
  src_baf <- file.path(rep_dir, "baf")
}

for (f in baf_files) {
  src_f <- file.path(src_baf, f)
  if (file.exists(src_f)) {
    file.copy(src_f, file.path(dir_ca_baf, f), overwrite = TRUE)
  } else {
    warning(paste("Missing CalicoST BAF file:", src_f))
  }
}

if (file.exists(meta_file)) {
  df_meta_ca <- read_parquet(meta_file)
  if ("index" %in% colnames(df_meta_ca)) {
    colnames(df_meta_ca)[colnames(df_meta_ca) == "index"] <- "barcode"
  } else {
    colnames(df_meta_ca)[1] <- "barcode"
  }
  x_vals <- if ("x" %in% colnames(df_meta_ca)) df_meta_ca$x else if ("x_coord" %in% colnames(df_meta_ca)) df_meta_ca$x_coord else df_meta_ca[[2]]
  y_vals <- if ("y" %in% colnames(df_meta_ca)) df_meta_ca$y else if ("y_coord" %in% colnames(df_meta_ca)) df_meta_ca$y_coord else df_meta_ca[[3]]
  df_tissue <- tibble(
    barcode = as.character(df_meta_ca$barcode),
    in_tissue = 1,
    x = as.numeric(x_vals),
    y = as.numeric(y_vals),
    pixel_row = as.numeric(y_vals),
    pixel_col = as.numeric(x_vals)
  )
  write.table(
    df_tissue,
    file.path(dir_ca_spatial, "tissue_positions_list.csv"),
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
} else {
  warning(paste("Missing CalicoST spatial metadata:", meta_file))
}

normalidx_file <- file.path(dir_ca_root, "normalidx_file.txt")
if (file.exists(normalidx_file)) {
  file.remove(normalidx_file)
}
message(">>> CalicoST inputs prepared in (ref-free mode): ", dir_ca_root)

# B. Ref-free ----
dir_cv <- file.path(out_root, "CopyVAE")
dir.create(dir_cv, recursive=TRUE, showWarnings=FALSE)
writeMM(mat_sparse, file.path(dir_cv, "counts.mtx"))
writeLines(rownames(mat_sparse), file.path(dir_cv, "genes.tsv"))
writeLines(colnames(mat_sparse), file.path(dir_cv, "barcodes.tsv"))

message("-> Done.")
