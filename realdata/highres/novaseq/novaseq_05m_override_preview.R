#!/usr/bin/env Rscript
# Preview of the IDC_Grade3 res0.3 compartment call after the manual override in
# novaseq_05_annotate.R (cluster 7 imm -> stroma), drawn from the current annotation
# table so it does not wait for the ~25 min annotation rerun.
#
#   spatial  three compartments, override applied
#   dotplot  12 clusters (labelled with the new call) x the 16 calling markers plus the
#            top genes (novaseq_05d top10) that decided the imm / stroma review
#
# Rationale: docs/IDC_Grade3_compartment_override.md
#
# Run (~2 min):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05m_override_preview.R

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(Seurat)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq", "processed_data", "clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")

OVERRIDE <- c(`7` = "stroma")          # same as novaseq_05_annotate.R
LEV <- c("epi", "imm", "stroma")
PAL <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
GENES <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  Immune     = c("CD3D", "CD3E", "JCHAIN", "MZB1", "MS4A2",
                 "PTPRC", "IGKC", "C1QA", "LYZ", "TPSB2"),
  Stromal    = c("COL1A1", "COL1A2", "LUM", "THBS2", "VWF", "PECAM1",
                 "POSTN", "MMP11", "APOD", "CXCL12", "CFD", "FABP4", "RGS5"))

## ---- new call per cluster ----
an <- fread(file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES)))
an[, cluster := as.character(cluster)]
an[, new_call := major_celltype]
an[cluster %in% names(OVERRIDE), new_call := OVERRIDE[cluster]]
print(an[, .(cluster, n_cells, old = major_celltype, new = new_call)][order(as.integer(cluster))])
stopifnot(identical(an[major_celltype == "epi", cluster], an[new_call == "epi", cluster]))

## ---- spatial ----
cells <- fread(file.path(CLUST, sprintf("%s_annotated_cells_res%s.csv", SAMPLE, RES)),
               select = c("x", "y", "cluster_frozen"))
cells[, compartment := factor(an$new_call[match(as.character(cluster_frozen), an$cluster)], levels = LEV)]
n <- cells[, .N, by = compartment][order(compartment)]
lab <- setNames(sprintf("%s (%s, %.1f%%)", n$compartment, format(n$N, big.mark = ","),
                        100 * n$N / sum(n$N)), n$compartment)
print(lab)
p_sp <- ggplot(cells, aes(x, y, colour = compartment)) +
  geom_point(size = 0.15, stroke = 0) +
  scale_colour_manual(values = PAL, labels = lab, name = "Compartment") +
  coord_fixed() + scale_y_reverse() + theme_void(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle(sprintf("%s res %s -- compartments after override (cluster 7 imm -> stroma)", SAMPLE, RES)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(OUT, sprintf("%s_res%s_override_spatial.png", SAMPLE, RES)), p_sp,
       width = 8, height = 8, dpi = 300)

## ---- dotplot: capped cells (2,000 per cluster, novaseq_05d) ----
b <- readRDS(file.path(OUT, sprintf("%s_res%s_markers_input.rds", SAMPLE, RES)))
g <- unlist(GENES, use.names = FALSE)
stopifnot(all(g %in% rownames(b$data)))
so <- CreateSeuratObject(CreateAssayObject(data = b$data[g, ]))
ord <- an[order(factor(new_call, levels = LEV), -n_cells)]
rl  <- setNames(sprintf("%s  %s%s", ord$cluster, ord$new_call,
                        ifelse(ord$cluster %in% names(OVERRIDE), " (was imm)", "")), ord$cluster)
so$row <- factor(unname(rl[b$idents[colnames(so)]]), levels = rev(unname(rl)))
# Seurat 5.1's grouped-feature DotPlot calls facet_grid(facets=), defunct in this ggplot2,
# so pass a plain vector and facet by gene group here.
p_dot <- DotPlot(so, features = g, group.by = "row", cols = c("grey90", "#B2182B"))
p_dot$data$panel <- factor(rep(names(GENES), lengths(GENES))[match(p_dot$data$features.plot, g)],
                           levels = names(GENES))
p_dot <- p_dot + facet_grid(cols = vars(panel), scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       title = sprintf("%s res %s -- clusters by new compartment", SAMPLE, RES),
       caption = "per group: first genes = calling markers, the rest = top genes from FindAllMarkers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8), plot.title = element_text(size = 9),
        plot.caption = element_text(size = 7, colour = "grey30"),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(OUT, sprintf("%s_res%s_override_dotplot.png", SAMPLE, RES)), p_dot,
       width = 11, height = 4.2, dpi = 300)
message("done -> ", OUT)
