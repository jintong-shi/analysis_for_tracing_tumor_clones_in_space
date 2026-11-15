#!/usr/bin/env Rscript
# =====================================================================
# Supplementary high-res figure (VisiumHD, AP5/pt10 prostate) — rebuilt from
# scratch using ONLY high_res_st/input_data. Every panel replotted from data;
# all derived values (histology collapse, meta-cell counts, deconv-vs-histology
# confusion, inferCNV obs downsample+clustering, subclone alignment) recomputed
# here. Sources ONLY utils/high_res_setup.R.
#
# Panels (letters = layout order; all titles removed as in the original):
#   a regular-Visium histology spatial      e deconv vs 8um histology heatmap
#   b HD 8um histology spatial              f per-unit QC violins (UMIs / genes)
#   c single-cell (RCTD) cell types         g res8 8um inferCNV heatmap (no group)
#   d n cells per meta-cell                 h epithelial(RCTD) inferCNV heatmap (no group)
#                                           i subclone HD-k2 over Visium-k3 (shared genome axis)
#                                           j 4-strategy pseudobulk inferCNV (was main Fig 5b)
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     paper/figures_final/high_res_st/scripts/supp_visiumHD.R
# =====================================================================
suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(grid)
}))

.this <- (function() { a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) normalizePath(sub("^--file=", "", m[1])) else normalizePath("supp_visiumHD.R") })()
SDIR <- dirname(.this); TDIR <- normalizePath(file.path(SDIR, "..")); UTILS <- normalizePath(file.path(TDIR, "..", "utils"))
IN   <- file.path(TDIR, "input_data", "visiumHD")
OUTF <- file.path(TDIR, "supp_figures"); dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)
OUTT <- file.path(TDIR, "data_tables");  dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS, "high_res_setup.R"))
wtab   <- function(df, nm) fwrite(as.data.frame(df), file.path(OUTT, nm))
wtab_gz <- function(df, nm) fwrite(as.data.frame(df), file.path(OUTT, nm))

FONT <- pick_font(); BASE <- 7; set.seed(42); PT_CAP <- 120000L
update_geom_defaults("text", list(family = FONT)); message("font: ", FONT)

## ------- palettes (verbatim from the frozen per-panel scripts) -------------
PAL_VIS <- c("Tumor"="#d80203","Benign"="#006cad","Lymphocytes/Inflammation"="#27b712",
             "Stroma"="#f2b8b1","Exclude"="#E0E0E0")
ORD_VIS <- c("Exclude","Stroma","Lymphocytes/Inflammation","Benign","Tumor")
PAL_HD  <- c("Tumor"="#d80203","Benign"="#006cad","Stroma"="#f2b8b1","Stroma-capsule"="#d47b70",
             "Stroma-intratumoural"="#e59a91","Stroma - blood vessel"="#f9dcd8","Stroma-nerve"="#a7709c",
             "Stroma - Lymphocytes"="#27b712","Stroma - TLS"="#5cc44a","Stoma - lymphocytes (benign)"="#8fd985")
ORD_HD  <- c("Stroma","Stroma-capsule","Stroma-intratumoural","Stroma - blood vessel","Stroma-nerve",
             "Stroma - Lymphocytes","Stroma - TLS","Stoma - lymphocytes (benign)","Benign","Tumor")
MAINTYPE_PAL <- c(Epithelial="#006cad", CAFs="#f2b8b1", SMCs="#e59a91", Endothelial="#f9dcd8",
                  PNS_glial="#a7709c", `B-cells`="#27b712", `T-cells`="#5cc44a", Myeloid="#1a8a0d",
                  Plasmablasts="#8fd985", Unknown="grey85")
SUB_COLS <- c(subclone_1="#b34d23", subclone_2="#9a76bf", subclone_3="#e77e71", subclone_4="#8172B3")
collapse_hist <- function(label) {
  l <- as.character(label); ll <- tolower(l); out <- l
  is_stroma <- grepl("^stroma|^stoma", ll)
  is_tum <- !is_stroma & grepl("tumour|tumor|cancer|idc", ll)
  is_ben <- !is_stroma & !is_tum & grepl("^benign", ll)
  out[is_tum] <- "Tumor"; out[is_ben] <- "Benign"; out
}
downsample <- function(df, cap = PT_CAP) if (nrow(df) > cap) df[sort(sample(nrow(df), cap)), , drop = FALSE] else df

theme_spatial <- function() theme_void(base_size = BASE) +
  theme(legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
        legend.key.size = unit(0.14, "cm"), legend.spacing.y = unit(0.5, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
        legend.position = "bottom", legend.direction = "horizontal",  # legend under each spatial panel -> full-width tissue (Jintong)
        plot.tag = element_text(size = 8, face = "bold"), plot.margin = margin(1, 1, 1, 1),
        plot.background = element_rect(fill = "white", color = NA))

# ============================================================ a: regular Visium
message("[a] regular-Visium histology")
coord_a <- fread(file.path(IN, "histology", "visium_tissue_positions.csv"))
annot_a <- fread(file.path(IN, "histology", "AP5_10_histology.csv"), check.names = FALSE)
da <- merge(coord_a, annot_a, by.x = "barcode", by.y = "Barcode")
da[, hist := collapse_hist(AP5_10_histology)]
da <- da[hist %in% names(PAL_VIS)]; da[, hist := factor(hist, levels = ORD_VIS)]; setorder(da, hist)
wtab(da[, .(barcode, x = pxl_col_in_fullres, y = pxl_row_in_fullres, hist)], "supp_visiumHD_panel_a_visium_histology.csv")
pa <- ggplot(da, aes(pxl_col_in_fullres, pxl_row_in_fullres, color = hist)) + geom_point(size = 0.42, stroke = 0) +
  scale_color_manual(values = PAL_VIS, name = "Histology (AP5_10)", drop = FALSE) +
  sq_scales(da$pxl_col_in_fullres, da$pxl_row_in_fullres) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.4))) + theme_spatial()  # 1 col: narrow, stays clear of b

# ================================================================= b: HD 8um
message("[b] HD 8um histology")
coord_b <- fread(file.path(IN, "histology", "hd8um_tissue_positions.csv"))
annot_b <- fread(file.path(IN, "histology", "AP5_10_HD_AL.csv"), check.names = FALSE)
db <- merge(coord_b, annot_b, by.x = "barcode", by.y = "Barcode")
db[, hist := collapse_hist(`Histology Annotations`)]
db <- db[hist %in% names(PAL_HD)]; db[, hist := factor(hist, levels = ORD_HD)]; setorder(db, hist)
wtab(db[, .(barcode, x = pxl_col_in_fullres, y = pxl_row_in_fullres, hist)], "supp_visiumHD_panel_b_hd8um_histology.csv")
db <- downsample(db)
pb <- ggplot(db, aes(pxl_col_in_fullres, pxl_row_in_fullres, color = hist)) + geom_point(size = 0.28, stroke = 0) +
  scale_color_manual(values = PAL_HD, name = "Histology Annotations", drop = FALSE) +
  sq_scales(db$pxl_col_in_fullres, db$pxl_row_in_fullres) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 1.4))) + theme_spatial()

# ====================================================== c: deconv cell types
message("[c] RCTD cell types")
dc <- fread(file.path(IN, "rctd_cells.csv"))
dc[, first_type := as.character(first_type)]; dc[is.na(first_type), first_type := "Unknown"]
lev_c <- c(setdiff(sort(unique(dc$first_type)), "Unknown"), "Unknown")
dc[, first_type := factor(first_type, levels = lev_c)]
setorder(dc, -first_type)                                    # Unknown underneath
wtab(dc[, .(x, y, first_type)], "supp_visiumHD_panel_c_rctd_celltypes.csv")
dc <- downsample(dc)
pal_c <- MAINTYPE_PAL[lev_c]; pal_c[is.na(pal_c)] <- "grey70"; names(pal_c) <- lev_c
pc <- ggplot(dc, aes(x, y, color = first_type)) + geom_point(size = 0.35, stroke = 0) +
  scale_color_manual(values = pal_c, name = "Cell type (RCTD)", drop = FALSE) +
  sq_scales(dc$x, dc$y) + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 1.4))) + theme_spatial()

# ====================================================== d: n cells per meta-cell
message("[d] n cells per meta-cell")
emb <- fread(file.path(IN, "epi_cell_embedding.csv"))            # cell, x, y
meta <- fread(file.path(IN, "meta_r40_Smax20.csv"))             # cell, pseudo_id
emb[, cell := as.character(cell)]; meta[, cell := as.character(cell)]
dd <- merge(emb, meta, by = "cell", all.x = TRUE)
dd <- dd[!is.na(pseudo_id) & pseudo_id != 0]
dd[, n_cells := .N, by = pseudo_id]
stopifnot(nrow(dd) > 0)
wtab(dd[, .(cell, x, y, pseudo_id, n_cells)], "supp_visiumHD_panel_d_metacell_ncells.csv")
dd <- downsample(dd)
pd <- ggplot(dd, aes(x, y, colour = n_cells)) + geom_point(size = 0.35, stroke = 0) +
  scale_colour_viridis_c(name = "n cells per meta-cell") + sq_scales(dd$x, dd$y) +
  guides(colour = guide_colourbar(barwidth = unit(1.3, "cm"), barheight = unit(0.24, "cm"))) + theme_spatial()

# ================================================ e: deconv vs histology heatmap
message("[e] deconv vs histology heatmap")
mc <- fread(file.path(IN, "deconv_vs_histology.csv"))          # Histology, Deconv (raw per-bin)
conf <- table(mc$Histology, mc$Deconv)
dfp  <- as.data.table(as.data.frame(conf)); setnames(dfp, c("Histology", "Deconv", "Freq"))
rs   <- dfp[, .(tot = sum(Freq)), by = Histology]
dfp  <- merge(dfp, rs, by = "Histology"); dfp[, prop := Freq / tot]
wtab(dfp[, .(Histology, Deconv, Freq, prop)], "supp_visiumHD_panel_e_deconv_vs_histology.csv")
pe <- ggplot(dfp, aes(Deconv, Histology, fill = prop)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Row-normalized\nproportion") +
  labs(x = "Deconvolution label", y = "Histology label (8um)") + theme_minimal(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 1), axis.text.y = element_text(size = BASE - 2),
        axis.title = element_text(size = BASE), legend.title = element_text(size = BASE),
        legend.text = element_text(size = BASE - 1), legend.key.size = unit(0.28, "cm"),
        panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA))

# ==================================================== f: per-unit QC (UMIs / genes)
# Ported from paper/figures/fig5_merge/panels/p2_qc.R. As the analysis unit coarsens
# (8um bin -> single cell -> meta-cell) per-unit depth rises: violin (width-scaled
# density) + boxplot, linear y, median annotated, n under each x label. EPI-only, so
# the three are like-for-like. Palette = fig5 strategy colours. Per-unit nCount/nFeature
# are extracted once by input_data/prep/prep_visiumHD_qc.R -> qc_per_unit.csv.
message("[f] per-unit QC (UMIs / genes)")
QC_COLS <- c("8µm bin" = "#C44E52", "Single cell" = "#DD8452", "Meta-cell" = "#55A868")
QC_LV   <- names(QC_COLS)
qc <- fread(file.path(IN, "qc_per_unit.csv"))
qc[, rep := factor(rep, levels = QC_LV)]
wtab(qc, "supp_visiumHD_panel_f_qc.csv")
summ <- qc[, .(n = .N, med_umi = round(median(nCount)), med_gene = round(median(nFeature))), by = rep]
setorder(summ, rep)
x_labs <- setNames(sprintf("%s\n(n=%s)", summ$rep, formatC(summ$n, format = "d", big.mark = ",")), summ$rep)
qc_panel <- function(yvar, ylab, med_col, ymax, ybreaks) {
  med <- data.frame(rep = factor(summ$rep, levels = QC_LV), y = summ[[med_col]])
  ggplot(qc, aes(rep, .data[[yvar]], fill = rep)) +
    geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "grey25", linewidth = 0.3) +
    geom_text(data = med, aes(rep, ymax * 0.96, label = formatC(y, format = "d", big.mark = ",")),
              inherit.aes = FALSE, hjust = 0.5, vjust = 1, size = 1.9) +  # median at top, no x-collision
    scale_fill_manual(values = QC_COLS, guide = "none") +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(breaks = ybreaks, labels = formatC(ybreaks, format = "d", big.mark = ",")) +
    coord_cartesian(ylim = c(0, ymax)) + labs(x = NULL, y = ylab) +
    theme_classic(base_size = BASE) +
    theme(axis.text.x  = element_text(size = BASE - 2.5, colour = "black", angle = 30, hjust = 1),
          axis.text.y  = element_text(size = BASE - 2, colour = "black"),
          axis.title.y = element_text(size = BASE - 1),
          plot.margin  = margin(2, 4, 2, 2), plot.background = element_rect(fill = "white", color = NA))
}
p_umi  <- qc_panel("nCount",   "UMIs per unit",  "med_umi",  15000, seq(0, 15000, 5000)) +
  labs(tag = "f") + theme(plot.tag = element_text(size = 8, face = "bold"))
p_gene <- qc_panel("nFeature", "Genes per unit", "med_gene",  7500, seq(0, 7500, 2500))
pqc <- wrap_plots(p_umi, p_gene, nrow = 1)

# ------------------------------------- inferCNV heatmaps (f/g/h shared chr axis)
INF <- file.path(IN, "infercnv"); SB <- file.path(IN, "subclone")
o8 <- read_infercnv(file.path(INF, "visiumHD_AP5_8um_run.final.infercnv_obj"))
oe <- read_infercnv(file.path(INF, "visiumHD_AP5_epithelial_fromRCTD_run.final.infercnv_obj"))
d2 <- readRDS(file.path(SB, "prep_r40_k2.rds"))
d3 <- readRDS(file.path(SB, "prep_visium_k3.rds"))
CHR_LV <- chr_levels_present(o8@gene_order[[1]], oe@gene_order[[1]], d2$chr_genes, drop_x = TRUE)  # autosomes only (Jintong)

chr_num <- ggplot2::as_labeller(function(x) sub("^chr", "", x))
theme_hm <- function(show_chr) theme_void(base_size = BASE) +
  theme(strip.text.x = if (show_chr) element_text(size = BASE - 2.5, margin = margin(b = 0.5)) else element_blank(),
        strip.text.y = element_blank(), panel.spacing.x = unit(0.3, "pt"), panel.spacing.y = unit(3, "pt"),
        panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.2),
        legend.title = element_text(size = BASE), legend.text = element_text(size = BASE - 1),
        legend.key.size = unit(0.24, "cm"), plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = margin(1, 1, 1, 1), plot.background = element_rect(fill = "white", color = NA))
hm_plot <- function(long, show_chr, show_legend, na_col = "grey50") ggplot(long, aes(gpos, row_id, fill = cnv)) + geom_raster() +
  facet_grid(cols = vars(chr), scales = "free_x", space = "fixed", drop = FALSE, labeller = chr_num) +
  cnv_fill_gradientn(na_col = na_col) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  guides(fill = if (show_legend) guide_colourbar(barwidth = unit(0.24, "cm"), barheight = unit(1.3, "cm")) else "none") +
  theme_hm(show_chr)

message("[g] res8 inferCNV"); long_f <- infercnv_long_obs(o8, CHR_LV)
message("[h] epithelial inferCNV"); long_g <- infercnv_long_obs(oe, CHR_LV)
wtab(as.data.table(long_f)[, .(chr, gpos, row_id, cnv = round(cnv, 4))], "supp_visiumHD_panel_g_res8_heatmap.csv.gz")
wtab(as.data.table(long_g)[, .(chr, gpos, row_id, cnv = round(cnv, 4))], "supp_visiumHD_panel_h_epithelial_heatmap.csv.gz")
pf <- hm_plot(long_f, show_chr = TRUE, show_legend = FALSE)
pg <- hm_plot(long_g, show_chr = FALSE, show_legend = FALSE)

# i: HD (k=2) over Visium (k=3); Visium subclone_2<->3 swapped to align to HD.
message("[i] subclone HD-k2 over Visium-k3")
DS_LV <- c("HD epi-meta-cells", "regular-Visium")
common <- intersect(colnames(d2$mat), colnames(d3$mat))
keep <- colnames(d2$mat) %in% common; genes <- colnames(d2$mat)[keep]
chr  <- factor(as.character(d2$chr_genes[keep]), levels = CHR_LV)
ord  <- order(chr); genes <- genes[ord]; chr <- chr[ord]
gpos <- ave(seq_along(chr), chr, FUN = seq_along)
swap23 <- function(s) { o <- s; o[s == "subclone_2"] <- "subclone_3"; o[s == "subclone_3"] <- "subclone_2"; o }
block <- function(mat, sub, dataset, cap, seed = 42L) {
  sub <- sub[rownames(mat)]
  set.seed(seed); if (nrow(mat) > cap) { i <- sort(sample(nrow(mat), cap)); mat <- mat[i, ]; sub <- sub[i] }
  o <- order(factor(sub, levels = paste0("subclone_", 1:4)))
  list(m = mat[o, genes, drop = FALSE], sub = sub[o], dataset = dataset) }
bhd  <- block(d2$mat, d2$sub,         DS_LV[1], 900L)
bvis <- block(d3$mat, swap23(d3$sub), DS_LV[2], 700L)
mk_long <- function(b) { ng <- ncol(b$m); nr <- nrow(b$m)
  data.frame(chr = rep(chr, nr), gpos = rep(gpos, nr), row_id = rep(seq_len(nr), each = ng),
             cnv = as.vector(t(b$m)), dataset = b$dataset) }
h_long  <- rbind(mk_long(bhd), mk_long(bvis))
h_strip <- rbind(data.frame(dataset = bhd$dataset,  row_id = seq_len(nrow(bhd$m)),  subclone = bhd$sub),
                 data.frame(dataset = bvis$dataset, row_id = seq_len(nrow(bvis$m)), subclone = bvis$sub))
h_long$dataset  <- factor(h_long$dataset,  levels = DS_LV); h_strip$dataset <- factor(h_strip$dataset, levels = DS_LV)
wtab(as.data.table(h_long)[, .(dataset, chr, gpos, row_id, cnv = round(cnv, 4))], "supp_visiumHD_panel_i_subclone_heatmap.csv.gz")
wtab(h_strip, "supp_visiumHD_panel_i_subclone_strip.csv")
message(sprintf("  common genes: %d | HD spots: %d | Visium spots: %d", length(common), nrow(bhd$m), nrow(bvis$m)))

hm_h <- ggplot(h_long, aes(gpos, row_id, fill = cnv)) + geom_raster() +
  facet_grid(rows = vars(dataset), cols = vars(chr), scales = "free", space = "free_y", drop = FALSE,
             labeller = labeller(chr = chr_num)) +
  cnv_fill_gradientn() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = unit(0.24, "cm"), barheight = unit(1.3, "cm"))) + theme_hm(FALSE)
strip_h <- ggplot(h_strip, aes(x = 1, y = row_id, fill = subclone)) + geom_raster() +
  facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = SUB_COLS, guide = "none") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) + labs(tag = "i") + theme_void(base_size = BASE) +
  theme(strip.text.y.left = element_text(angle = 90, size = BASE - 2, face = "bold"), strip.placement = "outside",
        panel.spacing.y = unit(3, "pt"), plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = margin(1, 1, 1, 18), plot.background = element_rect(fill = "white", color = NA))

# j: the 4 resolution strategies' pseudobulk inferCNV. Moved here from main Fig 5
# (Jintong 2026-10-02); the main figure keeps only the correlation bar that
# summarises it. Same chr axis convention as g/h/i -- autosomes only.
message("[j] 4-strategy pseudobulk inferCNV")
pbk    <- compute_pseudobulk(file.path(TDIR, "input_data"))
j_keep <- as.character(pbk$chr_split) %in% CHR_LV
j_chr  <- factor(as.character(pbk$chr_split)[j_keep], levels = CHR_LV)
j_mat  <- pbk$heat_mat[, j_keep, drop = FALSE]
j_ord  <- order(j_chr); j_chr <- j_chr[j_ord]; j_mat <- j_mat[, j_ord, drop = FALSE]
j_gpos <- ave(seq_along(j_chr), j_chr, FUN = seq_along)
j_long <- data.frame(chr = rep(j_chr, nrow(j_mat)), gpos = rep(j_gpos, nrow(j_mat)),
                     row_id = rep(seq_len(nrow(j_mat)), each = ncol(j_mat)),
                     cnv = as.vector(t(j_mat)))
wtab(data.table(strategy = LABEL_MAP[rownames(j_mat)][j_long$row_id], chr = j_long$chr,
                gpos = j_long$gpos, cnv = round(j_long$cnv, 4)),
     "supp_visiumHD_panel_j_pseudobulk_strategies.csv.gz")
message(sprintf("  strategies: %d | genes on autosomes: %d", nrow(j_mat), ncol(j_mat)))
pj <- hm_plot(j_long, show_chr = TRUE, show_legend = TRUE, na_col = "grey88") +
  scale_y_reverse(expand = c(0, 0), breaks = seq_len(nrow(j_mat)),
                  labels = LABEL_MAP[rownames(j_mat)]) +
  labs(tag = "j") +
  theme(axis.text.y = element_text(size = BASE - 1, hjust = 1, margin = margin(r = 1)),
        plot.margin = margin(1, 1, 1, 4))

## ---- standalone per-panel exports ----
save1 <- function(p, stem, w, h) save_panel_both(p, stem, OUTF, w, h, FONT)
save1(pa, "supp_visiumHD_panel_a_visium_histology", 3.2, 3.2)
save1(pb, "supp_visiumHD_panel_b_hd8um_histology",  3.2, 3.2)
save1(pc, "supp_visiumHD_panel_c_rctd_celltypes",   3.2, 3.2)
save1(pd, "supp_visiumHD_panel_d_metacell_ncells",  3.2, 3.2)
save1(pe, "supp_visiumHD_panel_e_deconv_vs_histology", 3.6, 3.0)
save1(pqc, "supp_visiumHD_panel_f_qc",              4.0, 2.4)
save1(pf, "supp_visiumHD_panel_g_res8_heatmap",     6.5, 1.4)
save1(pg, "supp_visiumHD_panel_h_epithelial_heatmap", 6.5, 1.4)
save1(hm_h, "supp_visiumHD_panel_i_subclone_heatmap", 6.5, 2.4)
save1(pj,   "supp_visiumHD_panel_j_pseudobulk_strategies", 6.5, 1.1)

## ---- compose (source-latest layout: a-d spatial w/ bottom legends + e|f QC row;
##      g/h/i shared-axis inferCNV block; H=10.5) ----
tag8 <- function(p, t) p + labs(tag = t)
lblank <- function(t) ggplot() + theme_void() + labs(tag = t) +
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.margin = margin(1, 1, 1, 1))
# a-d: legend under each spatial panel (theme_spatial legend.position="bottom") so the
# tissue uses full cell width. Row 1 = a|b, row 2 = c|d. c narrowed to 2 cols here so
# its legend fits under c and does not spill into d (Jintong).
ta <- tag8(pa, "a"); tb <- tag8(pb, "b"); td <- tag8(pd, "d")
tc <- tag8(pc, "c") + guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 1.4)))
abcd <- (ta | tb) / (tc | td) & theme(legend.title.position = "top")
# row 3 = e (deconv vs histology heatmap) | f (per-unit QC violins, pre-tagged "f"); f wider.
row_ef <- wrap_plots(tag8(pe, "e"), pqc, nrow = 1, widths = c(1, 1.5))
top <- wrap_plots(abcd, row_ef, ncol = 1, heights = c(5.0, 1.6))
# g/h/i in ONE flat 2-col grid so patchwork aligns the chr columns across rows
# (col 1 = panel letter / subclone strip; col 2 = the faceted heatmaps).
ghi <- wrap_plots(lblank("g"), pf, lblank("h"), pg, strip_h, hm_h, ncol = 2,
                  widths = c(0.06, 1), heights = c(1.0, 1.0, 2.2))
# shared Subclone legend under g/h/i, keyed to i's colour strip (fig5 g panel style, Jintong)
L_sub <- cowplot::get_legend(
  ggplot(data.frame(x = 1:3, s = factor(paste0("subclone_", 1:3)))) + geom_tile(aes(x, 1, fill = s)) +
    scale_fill_manual(name = "Subclone", values = SUB_COLS, labels = c("Clone 1", "Clone 2", "Clone 3")) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(size = BASE),
          legend.text = element_text(size = BASE - 1), legend.key.size = unit(0.24, "cm"), legend.margin = margin(0, 0, 0, 0)))
# free() releases the g/h/i block from cross-block left alignment (§11.1).
ghi_full <- wrap_plots(patchwork::free(ghi, type = "panel", side = "l"), wrap_elements(full = L_sub),
                       ncol = 1, heights = c(1, 0.05))
final <- wrap_plots(top, ghi_full, pj, ncol = 1, heights = c(6.6, 3.7, 0.95)) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

W <- 7.09; H <- 11.4
save_panel_both(final, "supp_high_res_visiumHD_compose", OUTF, W, H, FONT)
message("DONE supp_visiumHD -> ", OUTF)
