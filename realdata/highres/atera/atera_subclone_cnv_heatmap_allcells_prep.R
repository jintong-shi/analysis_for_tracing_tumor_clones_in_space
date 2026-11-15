#!/usr/bin/env Rscript
# All-epithelial subclone CNV heatmap - prep (step 1/2).
#
# Unlike atera_subclone_cnv_heatmap.R (6,001 representative cells from the
# subclone object), this pulls the FULL strict-QC inferCNV matrix and keeps
# EVERY epithelial (DCIS) cell present in it (84,247), each labelled by subclone
# (clean merged into noisy), plus the Benign cells (2,799) as a reference panel.
# Saves a compact RDS for the render step. Two-step split because infercnv and
# ComplexHeatmap live in different R envs (mirrors atera_bvc_prep/plot).
#
# Run with shared R (has infercnv):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/atera/code/atera_subclone_cnv_heatmap_allcells_prep.R

suppressWarnings(suppressMessages({
  library(infercnv)
  library(data.table)
}))

ATERA    <- "<DATA_ROOT>/high_res_ST/atera"
obj_path <- file.path(ATERA, "processed_data/run4_benign_ref_strict_qc/run.final.infercnv_obj")
hist_path<- file.path(ATERA, "data/histology_annotation.csv")
sc_path  <- file.path(ATERA, "processed_data/subclone/subclones.csv")
out_rds  <- file.path(ATERA, "processed_data/subclone/outputs",
                      "atera_subclone_cnv_input.rds")
SEED <- 42L

## Optional label surgery, same convention as the app's 22_/23_/24_ figure scripts:
##   --exclude subclone_1                 drop these subclones (their cells leave the panel)
##   --rename  subclone_2=A,subclone_3=C  rename; two old labels -> one new label MERGES them
args    <- commandArgs(trailingOnly = TRUE)
argval  <- function(f) if (f %in% args) args[[which(args == f) + 1L]] else ""
exclude <- if (nzchar(argval("--exclude"))) strsplit(argval("--exclude"), ",")[[1]] else character(0)
kv      <- if (nzchar(argval("--rename"))) strsplit(strsplit(argval("--rename"), ",")[[1]], "=") else list()
relabel <- setNames(vapply(kv, `[`, "", 2L), vapply(kv, `[`, "", 1L))

message("loading inferCNV object (~1.5G) ...")
o    <- readRDS(obj_path)
expr <- o@expr.data           # genes x cells
go   <- o@gene_order          # chr, start, stop
cells<- colnames(expr)

## ---- histology: Benign (ref) vs DCIS (epithelial obs) ----
hist <- fread(hist_path)
lab  <- ifelse(grepl("^dcis",   hist$histology), "DCIS",
        ifelse(grepl("^normal", hist$histology), "Benign", NA_character_))
cell2lab <- setNames(lab, hist$Cell.ID)
label <- cell2lab[cells]

## ---- subclone label per cell (expand saved_barcodes) ----
## The log is append-only, so a label can appear several times; the latest save wins.
sc <- fread(sc_path)
setorder(sc, timestamp)
sc <- sc[, .SD[.N], by = roi_label]
setorder(sc, roi_label)
bc2sub <- unlist(lapply(seq_len(nrow(sc)), function(i) {
  bcs <- strsplit(sc$saved_barcodes[i], ";", fixed = TRUE)[[1]]
  setNames(rep(sc$roi_label[i], length(bcs)), bcs)
}))
sub_all <- bc2sub[cells]

bad <- setdiff(c(exclude, names(relabel)), sc$roi_label)
if (length(bad)) stop("unknown subclone(s) in --exclude/--rename: ",
                      paste(bad, collapse = ", "), call. = FALSE)
drop_cell <- !is.na(sub_all) & sub_all %in% exclude
sub_all[sub_all %in% names(relabel)] <- relabel[sub_all[sub_all %in% names(relabel)]]

benign_idx <- which(label == "Benign")
dcis_idx   <- which(label == "DCIS" & !drop_cell)
if (length(exclude))
  message(sprintf("excluded %s: %d DCIS cells dropped",
                  paste(exclude, collapse = ","), sum(drop_cell & label == "DCIS")))
sub_dcis   <- sub_all[dcis_idx]
if (length(benign_idx) == 0L || length(dcis_idx) == 0L)
  stop("need both Benign and DCIS cells in the inferCNV matrix")
if (anyNA(sub_dcis))
  stop(sprintf("%d epithelial (DCIS) cells have no subclone label", sum(is.na(sub_dcis))))
message(sprintf("benign(ref)=%d  epithelial/DCIS=%d", length(benign_idx), length(dcis_idx)))
message("epithelial subclone breakdown:"); print(table(sub_dcis))

## ---- genome-ordered, chromosome-split columns ----
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
chr_levels <- chr_levels[chr_levels %in% unique(as.character(go$chr))]
chr_genes  <- factor(as.character(go$chr), levels = chr_levels)
gene_ord   <- order(chr_genes, go$start)
chr_genes  <- chr_genes[gene_ord]

build <- function(idx) {
  mm <- t(expr[gene_ord, idx, drop = FALSE])   # cells x genes
  rownames(mm) <- cells[idx]
  mm
}
mat_benign <- build(benign_idx)
mat_epi    <- build(dcis_idx)
sub_epi    <- sub_dcis

## random within-panel row order (subclone blocks re-formed at plot time via row_split)
set.seed(SEED)
mat_benign <- mat_benign[sample(nrow(mat_benign)), , drop = FALSE]
ord        <- sample(nrow(mat_epi))
mat_epi    <- mat_epi[ord, , drop = FALSE]
sub_epi    <- sub_epi[ord]

saveRDS(list(
  mat_benign = mat_benign,
  mat_epi    = mat_epi,
  sub_epi    = sub_epi,
  chr_genes  = chr_genes,
  n_benign   = nrow(mat_benign),
  n_epi      = nrow(mat_epi)
), out_rds)

## per-cell label table, for panels that need the labels without the 1 GB matrix
fwrite(data.table(cell_id = rownames(mat_epi), clone = sub_epi),
       sub("_cnv_input.rds$", "_cell_labels.csv", out_rds))
message("DONE prep -> ", out_rds)
