#!/usr/bin/env Rscript
# Deliverable 2b -- freeze one clustering resolution and annotate its clusters into the
# epithelial / immune / stromal compartments.
#
# Compartment calls come from novaseq_04_cluster.R's marker z-scores (canonical panel:
# epi EPCAM/KRT8/KRT18/KRT19/MUC1, imm CD3D/CD3E/JCHAIN/MZB1/MS4A2,
# stroma COL1A1/COL1A2/LUM/THBS2/VWF/PECAM1). Labels are written as `major_celltype`
# with the same epi/imm/stroma levels and palette as high_res_ST/atera, so the two
# platforms stay directly comparable.
#
# Clusters whose top-to-runner-up margin is below AMBIG_MARGIN are still given their
# best-scoring label, but are flagged in `compartment_ambiguous` so a mixed cluster
# cannot silently become part of the inferCNV reference.
#
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05_annotate.R --sample DCIS_IDC_Grade1 --res 0.1

suppressWarnings(suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork)
}))

args   <- commandArgs(trailingOnly = TRUE)
sample <- args[which(args == "--sample") + 1]
res    <- args[which(args == "--res") + 1]

PROJ <- "<DATA_ROOT>/high_res_ST/novaseq"
OUT  <- file.path(PROJ, "processed_data", "clustering", sample)

AMBIG_MARGIN <- 0.5
LEV     <- c("epi", "imm", "stroma")
PAL     <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
MAP     <- c(epithelial = "epi", immune = "imm", stromal = "stroma")
MARKERS <- c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1",
             "CD3D", "CD3E", "JCHAIN", "MZB1", "MS4A2",
             "COL1A1", "COL1A2", "LUM", "THBS2", "VWF", "PECAM1")

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

rc <- paste0("spatial_snn_res.", res)
log("reading ", sample, " object")
seu <- readRDS(file.path(OUT, paste0(sample, "_seurat.rds.gz")))
stopifnot(rc %in% colnames(seu@meta.data))

asg <- read.csv(file.path(OUT, paste0(sample, "_compartment_by_resolution.csv")))
asg <- asg[as.character(asg$resolution) == as.character(res), ]
asg$cluster <- as.character(asg$cluster)
asg$major_celltype <- unname(MAP[asg$compartment])
asg$ambiguous <- asg$margin < AMBIG_MARGIN

# Manual overrides after reviewing each cluster's FindAllMarkers top genes
# (results/strataMap_module_scores/<sample>_res<res>_top10_markers.csv, novaseq_05d).
# IDC_Grade3 res0.3 cluster 7: top genes APOD/C3/CXCL12/MMP2/CFD/CTSK (fibroblast /
# adipose stroma), PTPRC in 9% of cells; called imm only because its median 259 genes
# leaves every marker low and a 14% JCHAIN signal tops the z-score (2026-10-04).
OVERRIDE <- list(IDC_Grade3 = list(`0.3` = c(`7` = "stroma")))
ov <- OVERRIDE[[sample]][[as.character(res)]]
asg$override <- ""
if (length(ov)) {
  i <- match(names(ov), asg$cluster)
  asg$override[i] <- paste0(asg$major_celltype[i], "->", ov)
  asg$major_celltype[i] <- unname(ov)
}
stopifnot(setequal(asg$cluster, levels(seu[[rc]][, 1])))

cl <- as.character(seu[[rc]][, 1])
seu$cluster_frozen        <- factor(cl, levels = levels(seu[[rc]][, 1]))
seu$major_celltype        <- factor(asg$major_celltype[match(cl, asg$cluster)], levels = LEV)
seu$compartment_margin    <- asg$margin[match(cl, asg$cluster)]
seu$compartment_ambiguous <- asg$ambiguous[match(cl, asg$cluster)]
seu$annotated_resolution  <- res

log("compartment sizes:")
print(table(seu$major_celltype))
log("cells in an ambiguous cluster: ", sum(seu$compartment_ambiguous), " (",
    round(100 * mean(seu$compartment_ambiguous), 1), "%)")

## ---- tables ----
tbl <- asg[order(-asg$n_cells),
           c("cluster", "n_cells", "major_celltype", "compartment", "margin",
             "ambiguous", "override", "epithelial", "immune", "stromal")]
tbl$pct_cells <- round(100 * tbl$n_cells / sum(tbl$n_cells), 2)
write.csv(tbl, file.path(OUT, paste0(sample, "_annotation_res", res, ".csv")), row.names = FALSE)
print(tbl, row.names = FALSE)

md <- seu@meta.data
md$umap1 <- Embeddings(seu, "umap")[, 1]; md$umap2 <- Embeddings(seu, "umap")[, 2]
write.csv(md[, c("x", "y", "umap1", "umap2", "nCount_spatial", "nFeature_spatial",
                 "area_um2", "cluster_frozen", "major_celltype",
                 "compartment_margin", "compartment_ambiguous")],
          file.path(OUT, paste0(sample, "_annotated_cells_res", res, ".csv")))

saveRDS(seu, file.path(OUT, paste0(sample, "_seurat_annotated_res", res, ".rds.gz")),
        compress = "gzip")
log("annotated object saved")

## ---- figures ----
save_both <- function(p, base, w, h) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  log("DONE -> ", basename(base), ".{png,pdf}")
}
base <- file.path(OUT, paste0(sample, "_annot_res", res))

lab <- setNames(paste0(tbl$cluster, " (", tbl$major_celltype, ")"), tbl$cluster)
md$cluster_lab <- factor(lab[as.character(md$cluster_frozen)], levels = lab[order(names(lab))])

umap_ct <- ggplot(md, aes(umap1, umap2, colour = major_celltype)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = PAL) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_classic() + ggtitle(paste0(sample, " res ", res, " -- compartment"))
umap_cl <- ggplot(md, aes(umap1, umap2, colour = cluster_lab)) +
  geom_point(size = 0.15, stroke = 0) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  theme_classic() + labs(colour = "cluster") + ggtitle("clusters")
save_both(umap_ct | umap_cl, paste0(base, "_umap"), 14, 5.5)

sp_ct <- ggplot(md, aes(x, y, colour = major_celltype)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = PAL) +
  coord_fixed() + scale_y_reverse() + theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle(paste0(sample, " res ", res, " -- compartment"))
sp_cl <- ggplot(md, aes(x, y, colour = cluster_lab)) +
  geom_point(size = 0.15, stroke = 0) + coord_fixed() + scale_y_reverse() + theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(colour = "cluster") + ggtitle("clusters")
save_both(sp_ct | sp_cl, paste0(base, "_spatial"), 14, 9)

Idents(seu) <- "cluster_frozen"
ord <- tbl$cluster[order(factor(tbl$major_celltype, levels = LEV), -tbl$n_cells)]
Idents(seu) <- factor(as.character(Idents(seu)), levels = rev(ord))
p <- DotPlot(seu, features = intersect(MARKERS, rownames(seu))) + RotatedAxis() +
  ggtitle(paste0(sample, " res ", res, " -- clusters ordered by compartment")) +
  theme(axis.text = element_text(size = 8), plot.title = element_text(size = 9))
save_both(p, paste0(base, "_dotplot"), 10, 2 + 0.3 * length(ord))

log("ALL DONE -> ", OUT)
