#!/usr/bin/env Rscript
# =====================================================================
# Shared setup for the figures_final / high_res_st figures (Fig 5 main +
# supp VisiumHD + supp atera Xenium). This is the ONLY script the three
# high-res figure scripts are allowed to source.
#
# Provides:
#   - bundled Arial-metric font (utils/fonts/LiberationSans*) + pick_font()
#   - CNV colour helpers: cnv_colorramp() (ComplexHeatmap), cnv_fill_gradientn() (ggplot)
#   - CHR_ALL / chr_levels_present() / short_chr()      genome-axis ordering
#   - read_infercnv() / infercnv_pseudobulk_obs() / infercnv_long_obs()
#       recompute pseudobulk means and long heatmap frames FROM the raw
#       inferCNV objects (no frozen stat is read)
#   - spatial_lims() (pad-to-box) and sq_scales() (square coord) helpers
#   - shared palettes: CNV_BREAKS, CLONE, HD_COLS, LABEL_MAP
#
# Self-locating (works when source()'d); no project-path globals, no dir
# creation side effects.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
  # circlize is only needed by cnv_colorramp() (ComplexHeatmap panels, fig5_main.R);
  # loaded lazily so scripts that only use the ggplot path run without it.
})

## ---- locate self + bundled font ------------------------------------------
.high_res_dir <- (function() {
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile))
    return(dirname(normalizePath(fr$ofile, mustWork = TRUE)))
  d <- Sys.getenv("FIG_UTILS_DIR", ""); if (nzchar(d)) return(normalizePath(d))
  getwd()
})()
.hr_fonts <- file.path(.high_res_dir, "fonts")
if (dir.exists(.hr_fonts) &&
    length(list.files(.hr_fonts, pattern = "\\.ttf$", ignore.case = TRUE))) {
  if (Sys.getenv("XDG_DATA_HOME") == "")  Sys.setenv(XDG_DATA_HOME = .high_res_dir)
  if (Sys.getenv("XDG_CACHE_HOME") == "") Sys.setenv(XDG_CACHE_HOME = file.path(tempdir(), "fc_cache"))
  dir.create(Sys.getenv("XDG_CACHE_HOME"), showWarnings = FALSE, recursive = TRUE)
}

# first registered family in `prefer` (real Arial, then Liberation/Nimbus metric
# clones). If none is registered (e.g. systemfonts absent in stcnv-r-tuscan) fall
# back to "Helvetica" -> rendered via the base pdf() device as the base-14
# Arial-metric font (live text). Never silently uses DejaVu.
pick_font <- function(prefer = c("Arial", "Liberation Sans", "Nimbus Sans", "Nimbus Sans L")) {
  fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
  for (f in prefer) if (f %in% fams) return(f)
  "Helvetica"
}

# device helpers: for "Helvetica" use base pdf() (base-14, Arial-metric, no font
# file needed); for a real registered TTF use cairo_pdf (embeds it). PNG via ragg
# if present else cairo. draw_panel() draws either a ggplot or a grid grob.
open_pdf_dev <- function(path, w, h, font) {
  if (identical(font, "Helvetica")) grDevices::pdf(path, width = w, height = h)
  else grDevices::cairo_pdf(path, width = w, height = h)
}
open_png_dev <- function(path, w, h, res = 300) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(path, width = w, height = h, units = "in", res = res, background = "white")
  else grDevices::png(path, width = w, height = h, units = "in", res = res, bg = "white", type = "cairo")
}
draw_panel <- function(p) {
  if (inherits(p, c("grob", "gTree", "gList"))) grid::grid.draw(p) else print(p)
}
# save one panel as pdf+png; handles ggplot OR grid grob (ComplexHeatmap grabExpr).
save_panel_both <- function(p, stem, outdir, w, h, font) {
  open_pdf_dev(file.path(outdir, paste0(stem, ".pdf")), w, h, font); draw_panel(p); grDevices::dev.off()
  open_png_dev(file.path(outdir, paste0(stem, ".png")), w, h);        draw_panel(p); grDevices::dev.off()
}

## ---- CNV colour scale (identical bounds everywhere) ----------------------
CNV_BREAKS <- c(0.85, 1.0, 1.15)
CNV_COLS   <- c("#2166AC", "white", "#B2182B")
cnv_colorramp <- function(bounds = CNV_BREAKS) {
  if (!requireNamespace("circlize", quietly = TRUE)) stop("cnv_colorramp() needs the circlize package")
  circlize::colorRamp2(bounds, CNV_COLS)
}
# na_col marks genes the strategy did not recover; only the 4-strategy pseudobulk
# panel has any, and there it is the point (HD 8um shares 146 genes with the GT).
cnv_fill_gradientn <- function(name = "inferCNV", lim = c(0.85, 1.15), na_col = "grey50") {
  scale_fill_gradientn(colours = CNV_COLS, values = scales::rescale(CNV_BREAKS),
                       limits = lim, oob = scales::squish, breaks = CNV_BREAKS, name = name,
                       na.value = na_col)
}

## ---- genome-axis ordering ------------------------------------------------
CHR_ALL <- paste0("chr", c(1:22, "X", "Y"))
short_chr <- function(x) sub("^chr", "", x)           # "chr1" -> "1"
# union of chromosomes present across gene-chr vectors, in canonical order.
# drop_y=TRUE removes chrY (near-empty) as the supp shared axis does; drop_x=TRUE
# additionally removes chrX (autosomes-only axis, supp-visiumHD inferCNV heatmaps).
chr_levels_present <- function(..., drop_y = TRUE, drop_x = FALSE) {
  u <- unique(unlist(lapply(list(...), as.character)))
  lv <- if (drop_x) paste0("chr", 1:22) else if (drop_y) paste0("chr", c(1:22, "X")) else CHR_ALL
  lv[lv %in% u]
}

## ---- inferCNV object helpers (recompute from raw) ------------------------
read_infercnv <- function(path) {
  o <- readRDS(path)
  stopifnot(.hasSlot(o, "expr.data"), .hasSlot(o, "gene_order"),
            .hasSlot(o, "observation_grouped_cell_indices"))
  o
}
.obs_idx <- function(o) {
  idx <- as.integer(unlist(o@observation_grouped_cell_indices, use.names = FALSE))
  idx[!is.na(idx) & idx >= 1L & idx <= ncol(o@expr.data)]
}
# pseudobulk = mean inferCNV over the observation (tumour) cells, per gene.
infercnv_pseudobulk_obs <- function(o) {
  idx <- sort(unique(.obs_idx(o)))
  list(pb = rowMeans(o@expr.data[, idx, drop = FALSE]), go = o@gene_order)
}
# long ggplot frame: obs cells (downsampled to cap, ward.D2 clustered), genes
# ordered by chromosome then start; gpos = within-chr gene index.
infercnv_long_obs <- function(o, chr_lv, cap = 2500L, seed = 42L) {
  cnv <- o@expr.data; go <- o@gene_order
  chr <- factor(as.character(go[[1]]), levels = chr_lv)
  ord <- order(chr, go[[2]]); genes <- rownames(go)[ord]; chr <- chr[ord]
  nna <- !is.na(chr); genes <- genes[nna]; chr <- chr[nna]   # drop genes off the axis (e.g. chrX/Y) so no NA facet
  gpos <- ave(seq_along(chr), chr, FUN = seq_along)
  obs <- sort(unique(.obs_idx(o)))
  set.seed(seed); if (length(obs) > cap) obs <- sort(sample(obs, cap))
  m <- cnv[genes, obs, drop = FALSE]
  m <- m[, hclust(dist(t(m)), "ward.D2")$order, drop = FALSE]
  ng <- nrow(m); nc <- ncol(m)
  data.frame(chr = rep(chr, nc), gpos = rep(gpos, nc),
             row_id = rep(seq_len(nc), each = ng), cnv = as.vector(m))
}

## ---- spatial layout helpers ----------------------------------------------
# pad x/y to a common box aspect (w/h) so a scatter fills its cell WITHOUT
# distortion (returns padded x/y ranges). Used by the fig5 spatial panels.
spatial_lims <- function(x, y, box_asp) {
  xr <- range(x); yr <- range(y); xs <- diff(xr); ys <- diff(yr)
  if (xs / ys < box_asp) { e <- (ys * box_asp - xs) / 2; xr <- xr + c(-e, e) }
  else                   { e <- (xs / box_asp - ys) / 2; yr <- yr + c(-e, e) }
  list(x = xr, y = yr)
}
# square coord window (zoom, no censoring; ratio locked -> tissue never distorted),
# y reversed (image convention). Used by the supp VisiumHD spatial panels.
sq_scales <- function(x, y, pad = 1.02) {
  xr <- range(x); yr <- range(y); half <- max(diff(xr), diff(yr)) / 2 * pad
  xl <- c(mean(xr) - half, mean(xr) + half); yl <- c(mean(yr) - half, mean(yr) + half)
  list(scale_y_reverse(), coord_fixed(ratio = 1, xlim = xl, ylim = rev(yl), expand = FALSE))
}

## ---- shared palettes / label maps ----------------------------------------
CLONE   <- c("Clone 1" = "#b34d23", "Clone 2" = "#9a76bf", "Clone 3" = "#e77e71")
HD_COLS <- c("HD 8um (direct)" = "#C44E52", "HD epithelial (RCTD)" = "#DD8452",
             "HD pseudospot (r40_Smax20)" = "#55A868", "Normal Visium (GT)" = "#4C72B0")
LABEL_MAP <- c("HD 8um (direct)"            = "HD 8µm",
               "HD epithelial (RCTD)"       = "HD epithelial cells",
               "HD pseudospot (r40_Smax20)" = "HD epi-meta-cells",
               "Normal Visium (GT)"         = "regular-Visium")

## ---- pseudobulk across the 4 VisiumHD resolution strategies ---------------
# Was private to fig5_main.R. The strategy heatmap moved to the VisiumHD
# supplement while the correlation bar stayed in the main figure, so both
# scripts need it; cached per session because it reads 4 inferCNV objects.
pb_cache <- NULL
compute_pseudobulk <- function(IN) {
  if (!is.null(pb_cache)) return(pb_cache)
  ID <- file.path(IN, "visiumHD", "infercnv")
  OBJS <- list(
    "Normal Visium (GT)"         = file.path(ID, "visium_pt10_AP5_GT_run.final.infercnv_obj"),
    "HD 8um (direct)"            = file.path(ID, "visiumHD_AP5_8um_run.final.infercnv_obj"),
    "HD epithelial (RCTD)"       = file.path(ID, "visiumHD_AP5_epithelial_fromRCTD_run.final.infercnv_obj"),
    "HD pseudospot (r40_Smax20)" = file.path(ID, "visiumHD_AP5_pseudospot_r40_Smax20_run.final.infercnv_obj"))
  GT <- "Normal Visium (GT)"; HD <- setdiff(names(OBJS), GT)
  pbs <- lapply(OBJS, function(p) infercnv_pseudobulk_obs(read_infercnv(p)))
  gt_pb <- pbs[[GT]]$pb; gt_go <- pbs[[GT]]$go
  corr <- rbindlist(lapply(HD, function(l) {
    g <- intersect(names(gt_pb), names(pbs[[l]]$pb))
    data.table(strategy = l, n_shared = length(g), pearson = cor(gt_pb[g], pbs[[l]]$pb[g])) }))
  chr_col <- as.character(gt_go[[1]]); start_col <- gt_go[[2]]
  chr_lv <- CHR_ALL[CHR_ALL %in% unique(chr_col)]
  chrf <- factor(chr_col, levels = chr_lv); ord <- order(chrf, start_col)
  genes <- rownames(gt_go)[ord]; chr_split <- chrf[ord]
  mapbb <- function(pb) { out <- setNames(rep(NA_real_, length(genes)), genes)
    cm <- intersect(names(pb), genes); out[cm] <- pb[cm]; out }
  hm <- t(sapply(names(OBJS), function(l) mapbb(pbs[[l]]$pb)))
  rownames(hm) <- names(OBJS); colnames(hm) <- genes
  hm <- hm[c(HD, GT), , drop = FALSE]
  bnd <- as.numeric(quantile(abs(gt_pb - 1), 0.99)); cnv_bounds <- c(1 - bnd, 1, 1 + bnd)
  pb_cache <<- list(heat_mat = hm, chr_split = chr_split, cnv_bounds = cnv_bounds, corr = corr)
  pb_cache
}

## ---- tool concordance on subclone pseudobulk profiles ---------------------
# Reads <platform>_10mb_pseudobulk_profiles.csv (rows "<subclone>-<tool>", cols
# "<chr>:<bin>") produced by high_res_ST/code/subclone_tool_concordance.R, and
# returns the 3x3 tool-vs-tool Pearson r averaged over that platform's subclones.
# One matrix per platform: the subclones are the replicates, not separate panels.
CNV_TOOLS <- c("inferCNV", "CopyKAT", "CONICSmat")
tool_concordance <- function(IN, platform) {
  f <- file.path(IN, "tool_concordance", sprintf("%s_10mb_pseudobulk_profiles.csv", platform))
  d <- fread(f); P <- as.matrix(d[, -1L]); rownames(P) <- d[[1L]]
  subs <- sort(unique(sub("-[^-]*$", "", rownames(P))))
  rs <- lapply(subs, function(s) cor(t(P[paste(s, CNV_TOOLS, sep = "-"), , drop = FALSE])))
  M <- Reduce(`+`, rs) / length(rs)
  dimnames(M) <- list(CNV_TOOLS, CNV_TOOLS)
  list(r = M, n_sub = length(subs), n_bin = ncol(P))
}
