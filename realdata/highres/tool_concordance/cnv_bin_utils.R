#!/usr/bin/env Rscript
# cnv_bin_utils.R -- shared pieces for comparing inferCNV / CopyKAT / CONICSmat.
#
# The three tools emit different feature grids (inferCNV: genes with gene_order;
# CopyKAT: 220kb bins with chrom/chrompos; CONICSmat: genes with CONICS_axis), so
# everything is averaged onto a common fixed-width bin grid on chr1-22 before any
# comparison. Sourced by compare_cnvtools_similarity.R and
# subclone_tool_concordance.R.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

HIGHRES <- "<DATA_ROOT>/high_res_ST"
BIN_MB_DEFAULT <- 10
BIN_MB_SET     <- c(1, 5, 10)   # resolutions the analyses are run at
EPS            <- 1e-12         # matches task1_cnv_vs_dna 01_standardize_cnv_output_realdata.R

cnv_cfg <- local({
  cfg <- list(
    atera = list(
      base     = file.path(HIGHRES, "atera/processed_data/cnvtools/atera"),
      infercnv = file.path(HIGHRES, "atera/processed_data/run4_benign_ref_strict_qc/run.final.infercnv_obj")
    ),
    novaseq = list(
      base     = file.path(HIGHRES, "novaseq/processed_data/cnvtools/IDC_Grade3_ds90k"),
      infercnv = file.path(HIGHRES, "novaseq/processed_data/infercnv_IDC_Grade3_ds90k/infcnv/run.final.infercnv_obj")
    ),
    visiumhd = list(
      base     = file.path(HIGHRES, "visiumHD/processed_data/cnvtools/pseudo_r40_Smax20"),
      infercnv = file.path(HIGHRES, "visiumHD/processed_data/run_CNV/epi_pseudo_spot/r40_Smax20/formal_run/run.final.infercnv_obj")
    )
  )
  # atera/novaseq also exist as 20k/60k subsamples (CopyKAT + CONICSmat cap out at
  # hclust's 65,536 cells); those share the parent's inferCNV run, so only `base`
  # differs and the shared-cell intersect does the rest.
  for (n in c("20k", "60k")) {
    cfg[[paste0("atera_ds", n)]]   <- list(base = paste0(cfg$atera$base, "_ds", n),
                                           infercnv = cfg$atera$infercnv)
    cfg[[paste0("novaseq_ds", n)]] <- list(base = paste0(cfg$novaseq$base, "_ds", n),
                                           infercnv = cfg$novaseq$infercnv)
  }
  cfg
})

# chr as bare 1..22; a feature lands in bin <chr>:<floor(pos / bin_mb)>
bin_id <- function(chr, pos, bin_mb = BIN_MB_DEFAULT) {
  chr <- sub("^chr", "", as.character(chr))
  ok  <- chr %in% as.character(1:22)
  out <- rep(NA_character_, length(chr))
  out[ok] <- sprintf("%s:%d", chr[ok], floor(pos[ok] / (bin_mb * 1e6)))
  out
}

# mean of each cell's feature values within each bin -> bins x cells
to_bins <- function(mat, bins) {
  keep <- !is.na(bins)
  mat  <- mat[keep, , drop = FALSE]
  bins <- factor(bins[keep])
  # incidence matrix scaled to per-bin means, so one sparse multiply does it all
  inc <- sparseMatrix(i = as.integer(bins), j = seq_along(bins), x = 1,
                      dims = c(nlevels(bins), length(bins)))
  out <- as.matrix((inc %*% mat) / as.numeric(table(bins)))
  rownames(out) <- levels(bins)
  out
}

# CopyKAT passes its matrix through a data.frame, so make.names() rewrites any
# barcode R considers non-syntactic: atera's "aaaajgij-1" comes back as
# "aaaajgij.1", while visiumHD's "g1" and novaseq's "G3_353479" survive intact.
# Canonicalising every name set the same way lines the three tools back up.
canon <- function(x) {
  y <- make.names(x)
  if (anyDuplicated(y)) stop("make.names() collides barcodes; cannot match cells")
  y
}

# the three tools' matrices for one platform, each binned and canon-named
load_tool_bins <- function(platform, bin_mb = BIN_MB_DEFAULT) {
  p <- cnv_cfg[[platform]]
  if (is.null(p)) stop("platform must be one of: ", paste(names(cnv_cfg), collapse = ", "))

  # All three are put on the project's common log2 scale first, exactly as
  # realdata_unified_by_dataset/code/task1_cnv_vs_dna/scripts/
  # 01_standardize_cnv_output_realdata.R does it, and only then binned --
  # mean(log2(x)) != log2(mean(x)), so the transform has to come first.

  # inferCNV: expr.data is a ratio centred on 1 -> log2
  message("[", platform, " @", bin_mb, "Mb] inferCNV ...")
  o  <- readRDS(p$infercnv)
  ic <- to_bins(log2(o@expr.data + EPS), bin_id(o@gene_order$chr, o@gene_order$start, bin_mb))
  rm(o); gc(FALSE)

  # CopyKAT: read the saved fit object rather than any of the text dumps -- its
  # $CNAmat is the same matrix as _copykat_final_results_bin_by_cell.txt (values
  # verified identical) but without that file's truncated last line, and it is
  # already log2-scale, so it is taken as-is. Rows are 220kb bins with
  # chrom/chrompos; colnames are make.names()-mangled because CNAmat is a
  # data.frame, which canon() below undoes consistently across tools.
  message("[", platform, " @", bin_mb, "Mb] CopyKAT ...")
  fit <- readRDS(file.path(p$base, "COPYKAT/copykat_results.rds"))
  if (inherits(fit, "try-error") || is.null(fit$CNAmat))
    stop("CopyKAT run failed for ", platform, ": copykat_results.rds holds no CNAmat")
  ck <- to_bins(as.matrix(fit$CNAmat[, -(1:3)]),
                bin_id(fit$CNAmat$chrom, fit$CNAmat$chrompos, bin_mb))
  rm(fit); gc(FALSE)

  # CONICSmat: CONICS_track_matrix is already log2 (normMat uses log2), as-is
  message("[", platform, " @", bin_mb, "Mb] CONICSmat ...")
  cm_raw <- readRDS(file.path(p$base, "CONICSmat/CONICS_track_matrix.rds"))
  ax     <- fread(file.path(p$base, "CONICSmat/CONICS_axis.csv"))
  cm <- to_bins(cm_raw, bin_id(ax$chromosome_name, ax$start_position, bin_mb))
  rm(cm_raw, ax); gc(FALSE)

  colnames(ic) <- canon(colnames(ic))
  colnames(ck) <- canon(colnames(ck))
  colnames(cm) <- canon(colnames(cm))

  list(inferCNV = ic, CopyKAT = ck, CONICSmat = cm,
       keep = canon(scan(file.path(p$base, "input/keep.txt"), what = "", quiet = TRUE)))
}

# subclone labels Jintong defined from the inferCNV results, as barcode -> label
load_subclones <- function(platform) {
  base <- sub("_ds[0-9]+k$", "", platform)
  if (base == "atera") {
    d <- read.csv(file.path(HIGHRES, "atera/processed_data/subclone/subclones.csv"),
                  stringsAsFactors = FALSE)
    d <- d[d$mode == "subclone", ]
    d <- d[d$timestamp == max(d$timestamp), ]        # latest save wins
    d <- d[!duplicated(d$roi_label), ]
    bc <- strsplit(d$saved_barcodes, ";")
    labs <- rep(d$roi_label, lengths(bc))
    names(labs) <- unlist(bc)
  } else if (base == "novaseq") {
    d <- fread(file.path(HIGHRES, "novaseq/processed_data/infercnv_IDC_Grade3_ds90k",
                         "outputs/IDC_Grade3_ds90k_cnv_merged_m3.csv"))
    # Jintong merged subclone 3 into 2: the two are not distinct enough to keep
    # apart, leaving novaseq with two subclones.
    d[merged == 3L, merged := 2L]
    labs <- paste0("subclone_", d$merged); names(labs) <- d$barcode
  } else if (base == "visiumhd") {
    x <- readRDS(file.path(HIGHRES, "visiumHD/processed_data/run_CNV/epi_pseudo_spot",
                           "r40_Smax20/subclone_analysis/prep_r40_k2.rds"))
    labs <- as.character(x$sub); names(labs) <- x$obs_cells
  } else stop("no subclone definition wired for ", platform)

  if (anyDuplicated(names(labs))) stop("a barcode appears in more than one subclone")
  names(labs) <- canon(names(labs))
  labs
}
