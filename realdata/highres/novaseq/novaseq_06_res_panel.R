#!/usr/bin/env Rscript
# Legended UMAP + spatial panel for ONE resolution.
#
# The 10-panel grids from novaseq_04_cluster.R are drawn without legends (they would
# not fit), so the cluster colours there cannot be read off. This draws a single
# resolution at full size with an explicit legend giving, per cluster, its size and
# its marker-based compartment call.
#
# Reads only <sample>_cell_metadata.csv + <sample>_compartment_by_resolution.csv, so it
# does not need the multi-GB Seurat object and runs in under a minute.
#
# Run:
#   XDG_CACHE_HOME="$TMPDIR/fccache" \
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_06_res_panel.R --sample IDC_Grade3 --res 0.1

suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))

args   <- commandArgs(trailingOnly = TRUE)
sample <- args[which(args == "--sample") + 1]
res    <- args[which(args == "--res") + 1]

PROJ <- "<DATA_ROOT>/high_res_ST/novaseq"
OUT  <- file.path(PROJ, "processed_data", "clustering", sample)

# 31 distinguishable hues, enough for the largest resolution (res 1.0 -> 31 clusters)
PAL <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
         "#DC0000", "#7E6148", "#B09C85", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD",
         "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
         "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
         "#FB9A99", "#E31A1C", "#FDBF6F")
PAL_CT <- c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868")
MAP    <- c(epithelial = "epi", immune = "imm", stromal = "stroma")

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

rc <- paste0("spatial_snn_res.", res)
log("reading metadata for ", sample)
md <- read.csv(file.path(OUT, paste0(sample, "_cell_metadata.csv")), row.names = 1,
               check.names = FALSE)
stopifnot(rc %in% colnames(md))

asg <- read.csv(file.path(OUT, paste0(sample, "_compartment_by_resolution.csv")))
asg <- asg[as.character(asg$resolution) == as.character(res), ]
asg$cluster <- as.character(asg$cluster)
asg$ct <- unname(MAP[asg$compartment])
asg <- asg[order(as.integer(asg$cluster)), ]

n_tot <- nrow(md)
asg$label <- sprintf("%s - %s  (%s; %.1f%%)%s", asg$cluster, asg$ct,
                     format(asg$n_cells, big.mark = ","),
                     100 * asg$n_cells / n_tot,
                     ifelse(asg$margin < 0.5, " *", ""))

md$cluster <- factor(asg$label[match(as.character(md[[rc]]), asg$cluster)],
                     levels = asg$label)
md$compartment <- factor(asg$ct[match(as.character(md[[rc]]), asg$cluster)],
                         levels = c("epi", "imm", "stroma"))
pal <- setNames(PAL[seq_len(nrow(asg))], asg$label)

log(nrow(asg), " clusters at res ", res)

save_both <- function(p, base, w, h) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  log("DONE -> ", basename(base), ".{png,pdf}")
}

sub <- paste0("* = compartment margin < 0.5 (mixed cluster); ", n_tot, " cells")

u_cl <- ggplot(md, aes(umap_1, umap_2, colour = cluster)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = pal, name = "cluster") +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  theme_classic() + ggtitle(paste0(sample, " - res ", res, " - clusters"))
s_cl <- ggplot(md, aes(x, y, colour = cluster)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = pal, name = "cluster") +
  coord_fixed() + scale_y_reverse() + theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  ggtitle("spatial") + labs(caption = sub)
save_both(u_cl | s_cl, file.path(OUT, paste0(sample, "_res", res, "_clusters_legended")), 17, 8)

u_ct <- ggplot(md, aes(umap_1, umap_2, colour = compartment)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = PAL_CT) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_classic() + ggtitle(paste0(sample, " - res ", res, " - compartment"))
s_ct <- ggplot(md, aes(x, y, colour = compartment)) +
  geom_point(size = 0.15, stroke = 0) + scale_colour_manual(values = PAL_CT) +
  coord_fixed() + scale_y_reverse() + theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle("spatial") + labs(caption = sub)
save_both(u_ct | s_ct, file.path(OUT, paste0(sample, "_res", res, "_compartment_legended")), 15, 8)

write.csv(asg[, c("cluster", "ct", "compartment", "n_cells", "margin", "label")],
          file.path(OUT, paste0(sample, "_res", res, "_colour_key.csv")), row.names = FALSE)
key <- data.frame(cluster = asg$cluster, compartment = asg$ct, n_cells = asg$n_cells,
                  margin = asg$margin, colour = unname(pal))
print(key, row.names = FALSE)
log("ALL DONE -> ", OUT)
