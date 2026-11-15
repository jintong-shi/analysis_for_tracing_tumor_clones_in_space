#!/usr/bin/env Rscript
# subclone_tool_concordance.R -- how the three CNV tools render the subclones
# Jintong defined from the inferCNV results.
#
#   Rscript subclone_tool_concordance.R <platform> [bin_mb]   # default 10
#
# Subclone membership comes from the inferCNV-based definitions only; the tools
# are never re-clustered. Two profile figures, rows ordered subclone1-inferCNV,
# subclone1-CopyKAT, subclone1-CONICSmat, subclone2-... in both:
#
#   *_pseudo_subclone_heatmap.pdf  one row per (subclone, tool): the subclone's
#       cells averaged into a single pseudobulk CNV profile.
#   *_singlecell_subclone_heatmap.pdf  the same blocks, but one row per cell
#       (up to N_PER cells per subclone, the same cells reused across tools).
#
# All three tools arrive on the common log2 scale (see cnv_bin_utils.R), but
# their amplitudes still differ several-fold, so heatmap rows are z-scored across
# bins -- otherwise the figure mostly shows which tool is which.
#
# The genome-wide 3x3 tool similarity lives in compare_cnvtools_similarity.R.
# Here the tools are compared per subclone, at two levels that must not be
# confused with each other:
#
#   *_subclone_pseudobulk_similarity.pdf  3x3 panel per subclone, correlating the
#       subclone's PSEUDOBULK profiles -- do the tools agree on what the subclone is
#   *_within_subclone_percell_r.csv       the same tool pairs scored cell by cell
#       and then summarised by subclone -- do the tools agree on each single cell
#
# Averaging first removes the single-cell noise, so the pseudobulk numbers run
# far higher (0.88-0.94 vs 0.57-0.89 on atera). Reporting one as the other
# overstates or understates agreement by a wide margin.

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
  library(data.table)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "cnv_bin_utils.R"))

MIN_CELLS <- 10
N_PER     <- 60          # cells drawn per subclone for the single-cell figure
TOOLS     <- c("inferCNV", "CopyKAT", "CONICSmat")
OUTDIR    <- file.path(HIGHRES, "results", "subclone_tool_concordance")

args     <- commandArgs(trailingOnly = TRUE)
platform <- args[1]
bin_mb   <- if (length(args) > 1) as.numeric(args[2]) else BIN_MB_DEFAULT
tag      <- sprintf("%s_%gmb", platform, bin_mb)

labs <- load_subclones(platform)
M    <- load_tool_bins(platform, bin_mb)

cells <- Reduce(intersect, c(list(M$keep, names(labs)), lapply(M[TOOLS], colnames)))
bins  <- Reduce(intersect, lapply(M[TOOLS], rownames))
bins  <- bins[order(as.integer(sub(":.*", "", bins)), as.integer(sub(".*:", "", bins)))]

labs <- labs[cells]
tab  <- table(labs)
drop <- names(tab)[tab < MIN_CELLS]
if (length(drop)) {
  message("[", tag, "] dropping subclones with <", MIN_CELLS, " shared cells: ",
          paste(drop, collapse = ", "))
  cells <- cells[!labs %in% drop]; labs <- labs[cells]
}
subs <- sort(unique(labs))
message("[", tag, "] ", length(cells), " cells, ", length(subs), " subclones, ",
        length(bins), " bins")
print(table(labs))
if (length(subs) < 2) stop("need at least 2 subclones to compare")

mats <- lapply(M[TOOLS], function(m) m[bins, cells]); rm(M); gc(FALSE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
stub <- file.path(OUTDIR, tag)

pal      <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
ann_cols <- list(tool     = setNames(c("#4C72B0", "#DD8452", "#55A868"), TOOLS),
                 subclone = setNames(grDevices::hcl.colors(length(subs), "Set2"), subs))
chr      <- sub(":.*", "", bins)

# z-score each row, then clip so one bin cannot wash the map out
zclip <- function(X) {
  Z <- t(scale(t(X))); Z[!is.finite(Z)] <- 0
  lim <- quantile(abs(Z), 0.99, na.rm = TRUE)
  list(Z = pmax(pmin(Z, lim), -lim), lim = lim)
}

# ---- figure 1: pseudobulk, one row per (subclone, tool) --------------------
P <- do.call(rbind, lapply(subs, function(s) {
  cs <- cells[labs == s]
  t(sapply(TOOLS, function(t) rowMeans(mats[[t]][, cs, drop = FALSE])))
}))
rownames(P) <- as.vector(outer(TOOLS, subs, function(t, s) paste(s, t, sep = "-")))
fwrite(data.table(profile = rownames(P), P), paste0(stub, "_pseudobulk_profiles.csv"))

annP <- data.frame(tool     = factor(rep(TOOLS, length(subs)), levels = TOOLS),
                   subclone = factor(rep(subs, each = length(TOOLS)), levels = subs),
                   row.names = rownames(P))
zp <- zclip(P)
pheatmap(zp$Z, cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = seq_len(length(subs) - 1) * length(TOOLS),
         gaps_col = cumsum(rle(chr)$lengths),
         annotation_row = annP, annotation_colors = ann_cols,
         labels_col = ifelse(duplicated(chr), "", chr), fontsize_col = 6,
         breaks = seq(-zp$lim, zp$lim, length.out = 101), color = pal,
         main = sprintf("%s -- pseudo-subclone CNV profiles (%gMb bins, row z-scored)",
                        platform, bin_mb),
         filename = paste0(stub, "_pseudo_subclone_heatmap.pdf"),
         width = 14, height = 2 + 0.25 * nrow(P))

# ---- figure 2: single cells, same block order ------------------------------
set.seed(1234)
pick <- unlist(lapply(subs, function(s) {
  cs <- cells[labs == s]
  sort(sample(cs, min(N_PER, length(cs))))
}))
pick_sub <- labs[pick]

SC <- do.call(rbind, lapply(subs, function(s) {
  cs <- pick[pick_sub == s]
  do.call(rbind, lapply(TOOLS, function(t) t(mats[[t]][, cs, drop = FALSE])))
}))
blk <- unlist(lapply(subs, function(s) {
  n <- sum(pick_sub == s); rep(paste(s, TOOLS, sep = "-"), each = n)
}))
rownames(SC) <- make.unique(blk)

annS <- data.frame(tool     = factor(sub(".*-", "", blk), levels = TOOLS),
                   subclone = factor(sub("-[^-]*$", "", blk), levels = subs),
                   row.names = rownames(SC))
zs <- zclip(SC)
pheatmap(zs$Z, cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = cumsum(rle(blk)$lengths), gaps_col = cumsum(rle(chr)$lengths),
         annotation_row = annS, annotation_colors = ann_cols,
         show_rownames = FALSE,
         labels_col = ifelse(duplicated(chr), "", chr), fontsize_col = 6,
         breaks = seq(-zs$lim, zs$lim, length.out = 101), color = pal,
         main = sprintf("%s -- single-cell CNV profiles by subclone and tool (%gMb bins, row z-scored)",
                        platform, bin_mb),
         filename = paste0(stub, "_singlecell_subclone_heatmap.pdf"),
         width = 14, height = 3 + 0.012 * nrow(SC))

# ---- per-subclone breakdown of the per-cell rho ----------------------------
prs <- list(c("inferCNV", "CopyKAT"), c("inferCNV", "CONICSmat"), c("CopyKAT", "CONICSmat"))
colcor <- function(A, B) {
  a <- sweep(A, 2, colMeans(A)); b <- sweep(B, 2, colMeans(B))
  colSums(a * b) / (sqrt(colSums(a^2)) * sqrt(colSums(b^2)))
}
within <- rbindlist(lapply(prs, function(pr) {
  data.table(platform = platform, bin_mb = bin_mb, subclone = labs,
             pair = paste(pr, collapse = " vs "),
             r = colcor(mats[[pr[1]]], mats[[pr[2]]]))
}))[, .(n = .N, median_r = median(r, na.rm = TRUE)), by = .(platform, bin_mb, subclone, pair)]
fwrite(within, paste0(stub, "_within_subclone_percell_r.csv"))
print(dcast(within, subclone + n ~ pair, value.var = "median_r"))

# ---- 3x3 tool similarity, one panel per subclone ---------------------------
# This correlates the SUBCLONE PSEUDOBULK profiles -- the rows of P built above,
# one per (subclone, tool) -- and not the per-cell r summarised by subclone. The
# two answer different questions: averaging a subclone's cells first removes the
# single-cell noise, so this asks whether the tools agree on what the subclone
# IS, while the per-cell version asks whether they agree cell by cell. The
# per-cell numbers stay in *_within_subclone_percell_r.csv.
#
# Only same-subclone comparisons are shown -- subclone1-of-one-tool against
# subclone2-of-another mixes "do the tools agree" with "are the subclones
# different", which is a different question and not what this figure answers.
S <- rbindlist(lapply(subs, function(s) {
  M <- P[paste(s, TOOLS, sep = "-"), , drop = FALSE]
  rownames(M) <- TOOLS
  rbindlist(lapply(prs, function(pr) {
    r <- cor(M[pr[1], ], M[pr[2], ])
    data.table(subclone = s, a = c(pr[1], pr[2]), b = c(pr[2], pr[1]), r = r)
  }))[, rbind(.SD, data.table(a = TOOLS, b = TOOLS, r = 1), fill = TRUE)
      ][, subclone := s][]
}))
S[, `:=`(a = factor(a, levels = TOOLS), b = factor(b, levels = rev(TOOLS)))]
fwrite(dcast(S[a != b], subclone ~ paste(a, b, sep = "-"), value.var = "r"),
       paste0(stub, "_subclone_pseudobulk_tool_r.csv"))

ggplot(S, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.3f", r)), size = 2.6) +
  facet_wrap(~ subclone, nrow = 1) +
  # colour centred on 0.5, matching the genome-wide 3x3
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = .5, limits = c(0, 1), name = "Pearson r") +
  labs(x = NULL, y = NULL,
       title = sprintf("%s -- tool agreement on each subclone's pseudobulk profile", platform),
       subtitle = sprintf("Pearson r between subclone pseudobulk log2 profiles, %d %gMb bins",
                          length(bins), bin_mb)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(stub, "_subclone_pseudobulk_similarity.pdf"),
       width = 2 + 1.9 * length(subs), height = 3.4)
message("[", tag, "] done -> ", OUTDIR)
