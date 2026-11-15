#!/usr/bin/env Rscript

# Convert HATCHet2 per-clone segment output (best.seg.ucn + best.bbc.ucn) into
# Ding-style per-sample WES CNV truth CSVs for the metastatic_prostate_ffpe
# dataset.
#
# Two batches are written:
#   include_normal -> truth_links/wes_cnv/
#     T = u_normal*2 + sum_k u_clone_k * (A_k + B_k)
#     Reflects bulk WGS, diluted by stromal/normal admixture.
#   tumor_only     -> truth_links/wes_cnv_tumor_only/
#     T = sum_k (u_clone_k / sum_k u_clone_k) * (A_k + B_k)
#     Drops u_normal and renormalizes tumor clones. Removes purity dilution; if
#     no tumor clone is present (sum_k u_clone_k = 0) we set T=2 (neutral).
# log2_copy_ratio = log2(T / 2). CALL is +/-/0 by +-LOG2_THR.
# NUM_POINTS_COPY_RATIO = number of best.bbc.ucn bins fully within the segment.
#
# Output columns mirror Ding output: CONTIG, START, END, NUM_POINTS_COPY_RATIO,
# MEAN_LOG2_COPY_RATIO, CALL, sample. write.csv adds the quoted integer rowname
# column so the on-disk layout matches the existing Ding files exactly.

suppressPackageStartupMessages({
  library(data.table)
})

SEG_FILE     <- "<DATA_ROOT>/public_data/prostate_metastatic/pt10_FF/best.seg.ucn"
BBC_FILE     <- "<DATA_ROOT>/public_data/prostate_metastatic/pt10_FF/best.bbc.ucn"
BASE_DIR     <- "<DATA_ROOT>/realdata/datasets/metastatic_prostate_ffpe/truth_links"
LOG2_THR     <- 0.3

# Spatial-sample aliases: list of seq_ids whose segment-level total CN should
# be averaged. Single-element list = direct rename; multi-element = mean across
# replicates.
SPATIAL_ALIASES <- list(
  "pt10_A3" = c("NE712562", "NE780530"),
  "pt10_A4" = c("NE737576")
)

cn_total <- function(cn_string) {
  parts <- strsplit(cn_string, "|", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2L) return(NA_real_)
    a <- suppressWarnings(as.numeric(p[1]))
    b <- suppressWarnings(as.numeric(p[2]))
    a + b
  }, numeric(1))
}

derive_call <- function(log2, thr) {
  out <- rep("0", length(log2))
  out[is.na(log2)] <- NA_character_
  out[!is.na(log2) & log2 >=  thr] <- "+"
  out[!is.na(log2) & log2 <= -thr] <- "-"
  out
}

read_hatchet_table <- function(path) {
  hdr <- readLines(path, n = 1L)
  hdr <- sub("^#", "", hdr)
  cols <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
  dt <- fread(path, sep = "\t", skip = 1L, header = FALSE,
              col.names = cols, na.strings = c("", "NA"))
  dt
}

# Compute weighted-avg total CN under a given proportion-handling mode.
#   mode = "include_normal": use u_normal + u_clone* as-is (renormalize to 1).
#   mode = "tumor_only":     zero out u_normal, renormalize across tumor clones;
#                            if no tumor signal, fall back to T=2 (neutral).
compute_T_avg <- function(seg, mode) {
  cn_cols <- grep("^cn_", names(seg), value = TRUE)
  u_cols  <- sub("^cn_", "u_", cn_cols)
  stopifnot(all(u_cols %in% names(seg)))
  cn_tot_mat <- sapply(cn_cols, function(c) cn_total(seg[[c]]))
  if (is.null(dim(cn_tot_mat))) cn_tot_mat <- matrix(cn_tot_mat, ncol = length(cn_cols))
  u_mat <- as.matrix(seg[, ..u_cols])
  storage.mode(u_mat) <- "numeric"

  if (mode == "tumor_only") {
    norm_idx <- which(cn_cols == "cn_normal")
    if (length(norm_idx) == 1L) u_mat[, norm_idx] <- 0
  }
  row_sum <- rowSums(u_mat, na.rm = TRUE)
  zero_rows <- row_sum == 0
  row_sum[zero_rows] <- 1
  u_mat <- u_mat / row_sum
  T_avg <- rowSums(cn_tot_mat * u_mat, na.rm = TRUE)
  T_avg[zero_rows] <- 2  # no tumor signal -> neutral
  T_avg
}

# Attach NUM_POINTS_COPY_RATIO via interval overlap (bbc bins fully within seg).
attach_bin_counts <- function(seg, bbc) {
  bbc_iv <- bbc[, .(SAMPLE, CONTIG, bstart = START, bend = END)]
  seg_iv <- seg[, .(seg_id = .I, SAMPLE, CONTIG, sstart = START, send = END)]
  setkey(bbc_iv, SAMPLE, CONTIG, bstart, bend)
  setkey(seg_iv, SAMPLE, CONTIG, sstart, send)
  ov <- foverlaps(bbc_iv, seg_iv,
                  by.x = c("SAMPLE", "CONTIG", "bstart", "bend"),
                  by.y = c("SAMPLE", "CONTIG", "sstart", "send"),
                  type  = "within", nomatch = NULL)
  bin_counts <- ov[, .(NUM_POINTS_COPY_RATIO = .N), by = .(seg_id)]
  seg[, seg_id := .I]
  seg <- merge(seg, bin_counts, by = "seg_id", all.x = TRUE, sort = FALSE)
  seg[, seg_id := NULL]
  seg[is.na(NUM_POINTS_COPY_RATIO), NUM_POINTS_COPY_RATIO := 0L]
  seg
}

write_batch <- function(seg, out_dir, label) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ding_cols <- c("CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO",
                 "MEAN_LOG2_COPY_RATIO", "CALL", "sample")
  seq_ids <- sort(unique(seg$SAMPLE))
  cat(sprintf("\n[Batch=%s] writing %d per-seq_id files -> %s\n",
              label, length(seq_ids), out_dir))
  for (s in seq_ids) {
    sub <- seg[SAMPLE == s, .(CONTIG, START, END, NUM_POINTS_COPY_RATIO,
                              MEAN_LOG2_COPY_RATIO, CALL)]
    sub[, sample := s]
    setcolorder(sub, ding_cols)
    out_path <- file.path(out_dir, paste0(s, "_wes_cnv.csv"))
    write.csv(as.data.frame(sub), out_path)
    cat(sprintf("  [Save] %s | segs=%d (+%d -%d 0:%d)\n",
                basename(out_path), nrow(sub),
                sum(sub$CALL == "+", na.rm = TRUE),
                sum(sub$CALL == "-", na.rm = TRUE),
                sum(sub$CALL == "0", na.rm = TRUE)))
  }
  for (alias in names(SPATIAL_ALIASES)) {
    sids <- SPATIAL_ALIASES[[alias]]
    have <- sids[sids %in% seq_ids]
    miss <- setdiff(sids, have)
    if (length(miss) > 0L) {
      cat(sprintf("  [Warn] alias %s missing seq_ids: %s -- skipping\n",
                  alias, paste(miss, collapse = ",")))
      next
    }
    sub <- seg[SAMPLE %in% have,
               .(T_avg = mean(T_avg, na.rm = TRUE),
                 NUM_POINTS_COPY_RATIO = as.integer(round(mean(NUM_POINTS_COPY_RATIO, na.rm = TRUE)))),
               by = .(CONTIG, START, END)]
    sub[, MEAN_LOG2_COPY_RATIO := log2(T_avg / 2)]
    sub[is.infinite(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := NA_real_]
    sub[, CALL := derive_call(MEAN_LOG2_COPY_RATIO, LOG2_THR)]
    sub[, sample := alias]
    sub[, T_avg := NULL]
    setcolorder(sub, ding_cols)
    chr_order <- unique(seg$CONTIG)
    sub[, CONTIG := factor(CONTIG, levels = chr_order)]
    setorder(sub, CONTIG, START)
    sub[, CONTIG := as.character(CONTIG)]
    out_path <- file.path(out_dir, paste0(alias, "_wes_cnv.csv"))
    write.csv(as.data.frame(sub), out_path)
    cat(sprintf("  [Save] %s | replicates=%s | segs=%d (+%d -%d 0:%d)\n",
                basename(out_path), paste(have, collapse = "+"), nrow(sub),
                sum(sub$CALL == "+", na.rm = TRUE),
                sum(sub$CALL == "-", na.rm = TRUE),
                sum(sub$CALL == "0", na.rm = TRUE)))
  }
}

main <- function() {
  if (!file.exists(SEG_FILE)) stop("missing seg file: ", SEG_FILE)
  if (!file.exists(BBC_FILE)) stop("missing bbc file: ", BBC_FILE)

  seg_raw <- read_hatchet_table(SEG_FILE)
  bbc_raw <- read_hatchet_table(BBC_FILE)
  setnames(seg_raw, "CHR", "CONTIG")
  setnames(bbc_raw, "CHR", "CONTIG")

  # Bin counts depend only on segmentation, not on the CN mode -- compute once.
  seg_raw <- attach_bin_counts(seg_raw, bbc_raw)

  modes <- list(
    include_normal = file.path(BASE_DIR, "wes_cnv"),
    tumor_only     = file.path(BASE_DIR, "wes_cnv_tumor_only")
  )
  for (m in names(modes)) {
    seg_m <- copy(seg_raw)
    seg_m[, T_avg := compute_T_avg(seg_m, mode = m)]
    seg_m[, MEAN_LOG2_COPY_RATIO := log2(T_avg / 2)]
    seg_m[is.infinite(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := NA_real_]
    seg_m[, CALL := derive_call(MEAN_LOG2_COPY_RATIO, LOG2_THR)]
    write_batch(seg_m, modes[[m]], m)
  }
}

main()
