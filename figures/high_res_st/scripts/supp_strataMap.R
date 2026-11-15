#!/usr/bin/env Rscript
# =====================================================================
# Supplementary high-res figure (StrataMap, breast IDC grade 3) — the StrataMap
# counterpart of supp_atera.R, built from ONLY high_res_st/input_data/strataMap/.
# Every panel is replotted from those extracts; sources ONLY utils/high_res_setup.R.
#
# StrataMap is the platform the repo's folder names call "novaseq"; the figure uses
# the real name, so the panel titles here and Fig 5 f/g (labelled NovaSeq X) refer
# to the same data.
#
# Panels (letters = layout order):
# Row 1 is the clustering, row 2 is how the epithelium was picked out of it:
#   a nFeature spatial (YlOrRd, log10)
#   b res.0.3 cluster spatial
#   c res.0.3 cluster UMAP
#   d per-cell epithelial AddModuleScore (13 genes) by compartment call
#   e epithelial / immune / stromal compartment spatial
#   f EPCAM spatial (YlOrRd)
#   g per-cell CNV, the same 1500 cells in the same order in all three tools
#
# Row 2 was epithelium-vs-rest / EPCAM / per-cluster module scores until
# 2026-10-04, when Jintong swapped in d and e. Panel g reproduces (does not copy)
# results/cnvtool_similarity/percell_heatmap_novaseq_ds20k_10mb.pdf.
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     paper/figures_final/high_res_st/scripts/supp_strataMap.R
# =====================================================================
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2); library(patchwork) }))

.this <- (function() { a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) normalizePath(sub("^--file=", "", m[1])) else normalizePath("supp_strataMap.R") })()
SDIR <- dirname(.this); TDIR <- normalizePath(file.path(SDIR, "..")); UTILS <- normalizePath(file.path(TDIR, "..", "utils"))
IN   <- file.path(TDIR, "input_data", "strataMap")
OUTF <- file.path(TDIR, "supp_figures"); dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)
OUTT <- file.path(TDIR, "data_tables");  dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)
source(file.path(UTILS, "high_res_setup.R"))
wtab <- function(df, nm) fwrite(as.data.frame(df), file.path(OUTT, nm))

FONT <- pick_font(); base_size <- 6; PT_SP <- 0.07; PT_UM <- 0.07
update_geom_defaults("text", list(family = FONT)); message("font: ", FONT)

TOOLS <- c("inferCNV", "CopyKAT", "CONICSmat")

df <- fread(file.path(IN, "strataMap_supp_cells.csv"))
df[, cluster  := factor(cluster, levels = as.character(sort(as.integer(unique(cluster)))))]
df[, celltype := factor(celltype, levels = c("epi", "imm", "stroma"))]
# Image y grows downward, so plotting it upward puts the section on its head.
# Fig 5's novaseq panels negate y for the same reason.
df[, y := -y]
message("cells: ", nrow(df))

# atera's cluster palette, extended to the 12 res0.3 clusters
pal_res <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
             "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#FDB462", "#BC80BD")
names(pal_res) <- levels(df$cluster)
pal_cmp <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
ylord  <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")

# Three maps across one row, so each legend sits under its own panel rather than
# beside it: at a third of the page width a side legend would take most of the cell.
theme_sp <- theme_void(base_size = base_size, base_family = FONT) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(size = base_size, family = FONT),
        legend.text = element_text(size = base_size - 1, family = FONT), legend.key.size = unit(0.28, "cm"),
        plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(2, 2, 2, 2))
theme_um <- theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = base_size, family = FONT, margin = margin(t = 1)),
        axis.title.y = element_text(size = base_size, family = FONT, margin = margin(r = 1)),
        legend.title = element_text(size = base_size, family = FONT),
        legend.text = element_text(size = base_size - 1, family = FONT), legend.key.size = unit(0.28, "cm"),
        plot.background = element_rect(fill = "white", color = NA))

# Point layers stay vector, as in supp_atera. At 715,479 cells that is 4.3x atera's
# cell count, so the composed pdf runs to a few hundred MB; the png is the one to
# read on screen.

spatial_dim <- function(colvar, pal, legname)
  ggplot(df, aes(x = x, y = y, color = .data[[colvar]])) +
    geom_point(size = PT_SP, stroke = 0) + scale_color_manual(values = pal, name = legname, drop = FALSE) +
    coord_fixed() + guides(color = guide_legend(override.aes = list(size = 1.6), nrow = 2,
                                               title.position = "top")) + theme_sp
# Decade breaks only: the default picks five labels that collide on a 1.6 cm bar.
spatial_feat <- function(featvar, legname, trans = "identity") { d <- df[order(df[[featvar]])]
  brk <- if (trans == "log10") 10^(0:4)[10^(0:4) >= min(d[[featvar]]) & 10^(0:4) <= max(d[[featvar]])] else waiver()
  ggplot(d, aes(x = x, y = y, color = .data[[featvar]])) + geom_point(size = PT_SP, stroke = 0) +
    scale_color_gradientn(colours = ylord, trans = trans, name = legname, breaks = brk) + coord_fixed() +
    guides(color = guide_colorbar(barwidth = unit(1.6, "cm"), barheight = unit(0.25, "cm"),
                                  title.position = "top")) + theme_sp }
umap_dim <- function(colvar, pal, legname)
  ggplot(df, aes(x = umap1, y = umap2, color = .data[[colvar]])) +
    geom_point(size = PT_UM, stroke = 0) + scale_color_manual(values = pal, name = legname, drop = FALSE) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 1.6), ncol = 1)) + theme_um

p_nfeat_sp <- spatial_feat("nFeature", "nFeature", trans = if (min(df$nFeature) > 0) "log10" else "identity")
p_res_sp   <- spatial_dim("cluster", pal_res, "Cluster\n(res 0.3)")
p_res_umap <- umap_dim("cluster", pal_res, "Cluster\n(res 0.3)")
p_cmp_sp   <- spatial_dim("celltype", pal_cmp, "Compartment")
p_epcam_sp <- spatial_feat("EPCAM", "EPCAM")

## ---- d: epithelial module score per cell, by compartment call ---------------
# AddModuleScore with a 13-gene epithelial panel, independent of the 5 genes the
# call was made from (high_res_ST/novaseq/code/novaseq_05l_three_module_addmodulescore.R).
# Cells are the 2,000-per-cluster cap, so cluster 0 (210,750 cells) does not swamp it.
ams <- fread(file.path(IN, "strataMap_percell_three_module_ams.csv"))
ams[, frozen := factor(frozen, levels = names(pal_cmp))]
p_vio <- ggplot(ams, aes(frozen, epithelial, fill = frozen)) +
  geom_violin(linewidth = 0.15, scale = "width", colour = "grey25") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.15, fill = "white") +
  scale_fill_manual(values = pal_cmp, guide = "none") +
  labs(x = "Compartment", y = "Epithelial module score\n(AddModuleScore, 13 genes)") +
  theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.text = element_text(size = base_size - 1, colour = "black", family = FONT),
        axis.title = element_text(size = base_size, family = FONT),
        plot.background = element_rect(fill = "white", color = NA))

## ---- g: per-cell CNV, three tools, one shared cell order ---------------------
# The extract is already z-scored and already in the shared hclust order, so the
# panel only reshapes and draws it. Rows run top to bottom within each tool block.
cnv <- fread(file.path(IN, "strataMap_percell_cnv.csv"))
bin_cols <- setdiff(names(cnv), c("tool", "cell", "row_in_tool"))
bin_cols <- bin_cols[order(as.integer(sub(":.*", "", bin_cols)), as.integer(sub(".*:", "", bin_cols)))]
g_long <- melt(cnv, id.vars = c("tool", "row_in_tool"), measure.vars = bin_cols,
               variable.name = "bin", value.name = "z")
g_long[, chr := factor(sub(":.*", "", as.character(bin)),
                       levels = as.character(sort(unique(as.integer(sub(":.*", "", bin_cols))))))]
g_long[, gpos := as.integer(sub(".*:", "", as.character(bin)))]
g_long[, tool := factor(tool, levels = TOOLS)]
zlim <- as.numeric(quantile(abs(g_long$z), 0.99, na.rm = TRUE))
n_show <- max(cnv$row_in_tool)
message(sprintf("panel g: %d cells x %d bins x %d tools, z clipped at +/-%.1f",
                n_show, length(bin_cols), length(TOOLS), zlim))

p_cnv <- ggplot(g_long, aes(gpos, row_in_tool, fill = z)) + geom_raster() +
  facet_grid(rows = vars(tool), cols = vars(chr), scales = "free", space = "free_x", switch = "y") +
  scale_fill_gradientn(colours = c("#2166AC", "white", "#B2182B"), limits = c(-zlim, zlim),
                       oob = scales::squish, breaks = c(-zlim, 0, zlim),
                       labels = sprintf("%+.1f", c(-zlim, 0, zlim)), name = "z") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  guides(fill = guide_colorbar(barwidth = unit(0.25, "cm"), barheight = unit(1.4, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = base_size, base_family = FONT) +
  theme(strip.text.x = element_text(size = base_size - 2, margin = margin(b = 0.5)),
        strip.text.y.left = element_text(size = base_size, angle = 90, margin = margin(r = 1)),
        strip.placement = "outside",
        panel.spacing.x = unit(0.3, "pt"), panel.spacing.y = unit(2, "pt"),
        panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.2),
        legend.title = element_text(size = base_size), legend.text = element_text(size = base_size - 1),
        plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(1, 1, 1, 1))

## ---- data tables (plotted values) -------------------------------------------
wtab(df[, .(cell, x, y, nFeature)], "supp_strataMap_panel_a_nFeature_spatial.csv.gz")
wtab(df[, .(cell, x, y, cluster, celltype)], "supp_strataMap_panel_be_spatial_discrete.csv.gz")
wtab(df[, .(cell, umap1, umap2, cluster)], "supp_strataMap_panel_c_umap.csv.gz")
wtab(ams[, .(cell, compartment = frozen, epithelial_module = epithelial)], "supp_strataMap_panel_d_epi_module.csv.gz")
wtab(df[, .(cell, x, y, EPCAM)], "supp_strataMap_panel_f_EPCAM_spatial.csv.gz")
wtab(g_long[, .(tool, row_in_tool, chr, gpos, z = round(z, 4))], "supp_strataMap_panel_g_percell_cnv.csv.gz")

save_both <- function(p, stem, w, h) {
  grDevices::cairo_pdf(file.path(OUTF, paste0(stem, ".pdf")), width = w, height = h); print(p); dev.off()
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file.path(OUTF, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300, background = "white")
  else png(file.path(OUTF, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300, bg = "white", type = "cairo")
  print(p); dev.off(); message("wrote ", stem) }

asp  <- diff(range(df$y)) / diff(range(df$x))
sp_w <- 2.4; sp_h <- sp_w * asp + 0.5          # + room for the legend underneath

save_both(p_nfeat_sp, "supp_strataMap_panel_a_nFeature_spatial", sp_w, sp_h)
save_both(p_res_sp,   "supp_strataMap_panel_b_cluster_spatial",  sp_w, sp_h)
save_both(p_res_umap, "supp_strataMap_panel_c_cluster_umap",     2.9, 2.4)
save_both(p_vio,      "supp_strataMap_panel_d_epi_module_violin", 2.4, 2.4)
save_both(p_cmp_sp,   "supp_strataMap_panel_e_compartment_spatial", sp_w, sp_h)
save_both(p_epcam_sp, "supp_strataMap_panel_f_EPCAM_spatial",    sp_w, sp_h)
save_both(p_cnv,      "supp_strataMap_panel_g_percell_cnv",      7.09, 2.6)

## ---- composed figure --------------------------------------------------------
# Row 1 the clustering (two maps + its UMAP), row 2 how the epithelium was picked
# out of it (epithelial module score, compartment map, EPCAM map), row 3 the full-width per-cell CNV block.
design <- "
AAAAAAAABBBBBBBBCCCCCCCC
DDDDDDDDEEEEEEEEFFFFFFFF
GGGGGGGGGGGGGGGGGGGGGGGG
"
w_full <- 7.09
h_row <- (w_full / 3) * asp + 0.5
h_cnv <- 2.6
compose <- p_nfeat_sp + p_res_sp + p_res_umap + p_vio + p_cmp_sp + p_epcam_sp + p_cnv +
  plot_layout(design = design, heights = c(h_row, h_row, h_cnv)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
H <- 2 * h_row + h_cnv + 0.3
message(sprintf("compose: %.2f x %.2f in (%.0f x %.0f mm)", w_full, H, w_full * 25.4, H * 25.4))
save_both(compose, "supp_high_res_strataMap_compose", w_full, H)
message("DONE supp_strataMap -> ", OUTF)
