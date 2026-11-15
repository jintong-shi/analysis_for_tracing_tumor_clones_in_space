#!/usr/bin/env Rscript
# =====================================================================
# Supplementary high-res figure (atera Xenium breast) — rebuilt from scratch
# using ONLY high_res_st/input_data/atera/xenium_supp_cells.csv (the minimal
# per-cell extract of the Xenium object). Every panel is replotted from that
# table; sources ONLY utils/high_res_setup.R.
#
# Panels (letters = layout order):
#   a nFeature_Xenium spatial (YlOrRd, log10, rotated 90 CCW)
#   b res.0.1 spatial dimplot (rotated 90 CCW)
#   c res.0.1 UMAP
#   d major_celltype spatial dimplot (rotated 90 CCW)
#   e major_celltype UMAP
#   f EPCAM spatial (YlOrRd, rotated 90 CCW)
#   g EPCAM violin by major cell type
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/scripts/supp_atera.R
# =====================================================================
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2); library(patchwork) }))

.this <- (function() { a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) normalizePath(sub("^--file=", "", m[1])) else normalizePath("supp_atera.R") })()
SDIR <- dirname(.this); TDIR <- normalizePath(file.path(SDIR, "..")); UTILS <- normalizePath(file.path(TDIR, "..", "utils"))
IN <- file.path(TDIR, "input_data", "atera")
OUTF <- file.path(TDIR, "supp_figures"); dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)
OUTT <- file.path(TDIR, "data_tables");  dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS, "high_res_setup.R"))
wtab <- function(df, nm) fwrite(as.data.frame(df), file.path(OUTT, nm))

FONT <- pick_font(); base_size <- 6; PT_SP <- 0.3; PT_UM <- 0.15
update_geom_defaults("text", list(family = FONT)); message("font: ", FONT)

df <- fread(file.path(IN, "xenium_supp_cells.csv"))
df[, res01 := factor(as.character(res01), levels = as.character(0:10))]
df[, celltype := factor(celltype, levels = c("epi", "imm", "stroma"))]
message("cells: ", nrow(df))

pal_res <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
             "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#FDB462"); names(pal_res) <- as.character(0:10)
pal_ct  <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
ylord   <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")

theme_sp <- theme_void(base_size = base_size, base_family = FONT) +
  theme(legend.position = "right", legend.title = element_text(size = base_size, family = FONT),
        legend.text = element_text(size = base_size - 1, family = FONT), legend.key.size = unit(0.28, "cm"),
        plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(2, 2, 2, 2))
theme_um <- theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = base_size, family = FONT, margin = margin(t = 1)),
        axis.title.y = element_text(size = base_size, family = FONT, margin = margin(r = 1)),
        legend.title = element_text(size = base_size, family = FONT),
        legend.text = element_text(size = base_size - 1, family = FONT), legend.key.size = unit(0.28, "cm"),
        plot.background = element_rect(fill = "white", color = NA))

# rotated 90 CCW = aes(x=y, y=x) with un-reversed axes; tissue aspect preserved
spatial_dim <- function(colvar, pal, legname) ggplot(df, aes(x = y, y = x, color = .data[[colvar]])) +
  geom_point(size = PT_SP, stroke = 0) + scale_color_manual(values = pal, name = legname, drop = FALSE) +
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 2))) + theme_sp
spatial_feat <- function(featvar, legname, trans = "identity") { d <- df[order(df[[featvar]])]
  ggplot(d, aes(x = y, y = x, color = .data[[featvar]])) + geom_point(size = PT_SP, stroke = 0) +
    scale_color_gradientn(colours = ylord, trans = trans, name = legname) + coord_fixed() +
    guides(color = guide_colorbar(barwidth = unit(0.25, "cm"), barheight = unit(1.6, "cm"))) + theme_sp }
umap_dim <- function(colvar, pal, legname) ggplot(df, aes(x = umap1, y = umap2, color = .data[[colvar]])) +
  geom_point(size = PT_UM, stroke = 0) + scale_color_manual(values = pal, name = legname, drop = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") + guides(color = guide_legend(override.aes = list(size = 2))) + theme_um

p_nfeat_sp <- spatial_feat("nFeature", "nFeature", trans = if (min(df$nFeature) > 0) "log10" else "identity")
p_epcam_sp <- spatial_feat("EPCAM", "EPCAM")
p_res_umap <- umap_dim("res01", pal_res, "Cluster")
p_res_sp   <- spatial_dim("res01", pal_res, "Cluster")
p_ct_umap  <- umap_dim("celltype", pal_ct, "Cell type")
p_ct_sp    <- spatial_dim("celltype", pal_ct, "Cell type")
p_vln <- ggplot(df, aes(x = celltype, y = EPCAM, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) + scale_fill_manual(values = pal_ct, guide = "none") +
  labs(x = NULL, y = "EPCAM (log-normalized)") + theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.text = element_text(size = base_size - 1, color = "black", family = FONT),
        axis.title.y = element_text(size = base_size, family = FONT, margin = margin(r = 1)),
        plot.background = element_rect(fill = "white", color = NA))

## data tables (plotted values) ----
wtab(df[, .(cell, plot_x = y, plot_y = x, nFeature, EPCAM)], "supp_atera_panel_af_spatial_continuous.csv")
wtab(df[, .(cell, plot_x = y, plot_y = x, res01, celltype)], "supp_atera_panel_bd_spatial_discrete.csv")
wtab(df[, .(cell, umap1, umap2, res01, celltype)], "supp_atera_panel_ce_umap.csv")
wtab(df[, .(cell, celltype, EPCAM)], "supp_atera_panel_g_epcam_violin.csv")

save_both <- function(p, stem, w, h) {
  grDevices::cairo_pdf(file.path(OUTF, paste0(stem, ".pdf")), width = w, height = h); print(p); dev.off()
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file.path(OUTF, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300, background = "white")
  else png(file.path(OUTF, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300, bg = "white", type = "cairo")
  print(p); dev.off(); message("wrote ", stem) }

xr <- diff(range(df$x)); yr <- diff(range(df$y)); asp <- xr / yr
sp_w <- 4.1; sp_h <- (sp_w - 0.9) * asp + 0.3

save_both(p_nfeat_sp, "supp_atera_panel_a_nFeature_spatial", sp_w, sp_h)
save_both(p_res_sp,   "supp_atera_panel_b_res01_spatial",    sp_w, sp_h)
save_both(p_res_umap, "supp_atera_panel_c_res01_umap",       2.9, 2.4)
save_both(p_ct_sp,    "supp_atera_panel_d_celltype_spatial", sp_w, sp_h)
save_both(p_ct_umap,  "supp_atera_panel_e_celltype_umap",    2.9, 2.4)
save_both(p_epcam_sp, "supp_atera_panel_f_EPCAM_spatial",    sp_w, sp_h)
save_both(p_vln,      "supp_atera_panel_g_EPCAM_violin",     2.4, 2.4)

## composed figure (24-col grid: maps left, summaries right; top-right empty) ----
design <- "
AAAAAAAAAAAAAA##########
BBBBBBBBBBBBBB#CCCCCCCCC
DDDDDDDDDDDDDD#EEEEEEEEE
FFFFFFFFFFFFFF#GGGGGGGGG
"
w_full <- 7.09; leg_w <- 0.72; h_row <- ((14 / 24) * w_full - leg_w) * asp
compose <- p_nfeat_sp + p_res_sp + p_res_umap + p_ct_sp + p_ct_umap + p_epcam_sp + p_vln +
  plot_layout(design = design, heights = rep(h_row, 4)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
save_both(compose, "supp_high_res_atera_compose", w_full, 4 * h_row + 0.3)
message("DONE supp_atera -> ", OUTF)
