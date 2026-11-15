#!/usr/bin/env Rscript
# Do the classic epithelial markers actually separate the res0.3 epi compartment from
# the other two? Only epithelial markers are plotted -- the question is whether epi
# scores high and imm/stroma score low, not whether each compartment has its own genes.
#
# Reads the capped matrix novaseq_05d step 'extract' wrote (2,000 cells per res0.3
# cluster, log-normalised), so the 7 GB object is not opened again. The cap matters
# here: cluster 0 has 210,750 cells and cluster 11 has 641, and an unweighted mean over
# the raw compartment would be cluster 0's mean with a rounding error.
#
# Colour is the group mean divided by the largest group mean FOR THAT GENE, not a
# z-score. Across three groups a z-score degenerates towards +-0.707/+-1.22 and the
# colour becomes little more than the sign, which is what produced the CPB1 artefact
# earlier in this analysis. Relative-to-max keeps a real magnitude: 1.0 is the highest
# of the groups shown and 0.1 means a tenth of it.
#
#   a  markers x the three frozen compartments
#   b  markers x the twelve clusters -- a compartment that is a mixture shows up here
#      and not in a
#   c  per-cell mean of the core panel, by cluster
#
# Run:
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05k_epi_markers_by_compartment.R

suppressWarnings(suppressMessages({
  library(data.table); library(Matrix); library(ggplot2); library(patchwork)
}))

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"; RES <- "0.3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# classic epithelial markers, grouped by what they are evidence of. The five the
# compartment call itself used are marked in the panel title so the circularity is
# visible: EPCAM, KRT8, KRT18, KRT19, MUC1.
PANEL <- list(
  "Pan-epithelial"  = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7"),
  "Junction"        = c("CLDN3", "CLDN4", "CLDN7", "ELF3"),
  "Luminal breast"  = c("ESR1", "GATA3", "FOXA1", "TFAP2A", "MUC1"),
  "Basal / myoepi"  = c("KRT5", "KRT14", "KRT17", "TP63", "MIR205HG"))
USED <- c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1")
CORE <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN4", "ELF3")

b   <- readRDS(file.path(OUT, sprintf("%s_res%s_markers_input.rds", SAMPLE, RES)))
ann <- fread(file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES)))
ann[, cluster := as.character(cluster)]
cl_lv <- as.character(sort(as.integer(ann$cluster)))
cl    <- factor(b$idents, levels = cl_lv)
froz  <- setNames(ann$major_celltype[match(cl_lv, ann$cluster)], cl_lv)

genes <- unlist(PANEL, use.names = FALSE)
miss  <- setdiff(genes, rownames(b$data))
if (length(miss)) log("not in the matrix, dropped: ", paste(miss, collapse = " "))
genes <- intersect(genes, rownames(b$data))
grp_of <- unlist(lapply(names(PANEL), function(k)
  setNames(rep(k, length(PANEL[[k]])), PANEL[[k]])))
X <- b$data[genes, , drop = FALSE]
log("cells: ", ncol(X), "  genes: ", nrow(X), "  clusters: ", nlevels(cl))

# mean expression and % expressing, for an arbitrary grouping of the cells
tab <- function(g) {
  g <- factor(g)
  m <- sapply(levels(g), function(k) Matrix::rowMeans(X[, g == k, drop = FALSE]))
  p <- sapply(levels(g), function(k) Matrix::rowMeans(X[, g == k, drop = FALSE] > 0)) * 100
  d <- data.table(gene = rep(rownames(m), ncol(m)),
                  group = rep(colnames(m), each = nrow(m)),
                  mean = as.vector(m), pct = as.vector(p))
  d[, rel := mean / max(mean), by = gene]        # 1.0 = the highest group for this gene
  d[, panel := factor(grp_of[gene], levels = names(PANEL))]
  d[, gene := factor(gene, levels = genes)]
  d[]
}

CMP_LV <- c("epi", "imm", "stroma")
t_cmp <- tab(froz[as.character(cl)])
t_cmp[, group := factor(group, levels = CMP_LV)]
t_cl  <- tab(as.character(cl))
t_cl[, group := factor(group, levels = cl_lv)]

fwrite(t_cmp, file.path(OUT, sprintf("%s_res%s_epi_markers_by_compartment.csv", SAMPLE, RES)))

cat("\n=== core pan-epithelial panel, by frozen compartment ===\n")
print(dcast(t_cmp[gene %in% CORE], gene ~ group, value.var = "pct")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 1) else z)])
cat("(% of cells expressing)\n")

th <- theme_bw(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(linewidth = .15, colour = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6.5),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 6.5, face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        legend.key.size = unit(0.26, "cm"), legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7))

# the five genes the compartment call was made from are starred, so a reader can see
# which columns are circular
star <- function(g) ifelse(as.character(g) %in% USED, paste0(as.character(g), " *"),
                           as.character(g))

dot <- function(d, ttl, sub, ylab) {
  d <- copy(d)[, glab := factor(star(gene), levels = star(levels(gene)))]
  ggplot(d, aes(glab, group, size = pct, colour = rel)) +
    geom_point(stroke = 0) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    scale_size_continuous(range = c(0.1, 3.4), limits = c(0, 100),
                          breaks = c(10, 50, 90), name = "% expressing") +
    scale_colour_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B"),
                           limits = c(0, 1), name = "mean / max") +
    labs(title = ttl, subtitle = sub, x = NULL, y = ylab) + th
}

p_a <- dot(t_cmp, "a  Epithelial markers by res0.3 frozen compartment",
           "* = one of the five genes the compartment call itself used", "frozen call")
lab_cl <- sprintf("%s  %s (%s)", cl_lv, froz[cl_lv],
                  formatC(ann$n_cells[match(cl_lv, ann$cluster)], big.mark = ",", format = "d"))
t_cl[, group := factor(lab_cl[match(as.character(group), cl_lv)], levels = rev(lab_cl))]
p_b <- dot(t_cl, "b  The same, by cluster",
           "a compartment that is a mixture of lineages shows it here and not in a",
           "cluster  call (n)")

## ---- c: per-cell mean of the core panel --------------------------------------
sc <- Matrix::colMeans(b$data[intersect(CORE, rownames(b$data)), , drop = FALSE])
s  <- data.table(cluster = as.character(cl), froz = froz[as.character(cl)], score = sc)
s[, lab := factor(lab_cl[match(cluster, cl_lv)], levels = rev(lab_cl))]
s[, froz := factor(froz, levels = CMP_LV)]
med <- s[froz == "epi", median(score)]
p_c <- ggplot(s, aes(lab, score, fill = froz)) +
  geom_hline(yintercept = med, linetype = 2, colour = "grey40", linewidth = .35) +
  geom_boxplot(outlier.shape = NA, linewidth = .25, width = .66) +
  scale_fill_manual(values = c(epi = "#C44E52", imm = "#4C72B0", stroma = "#55A868"),
                    name = "frozen call") +
  coord_flip() +
  labs(title = "c  Per-cell mean of the core pan-epithelial panel",
       subtitle = sprintf("%s; dashed = median of the epi compartment",
                          paste(CORE, collapse = " ")),
       x = NULL, y = "mean log-normalised expression") + th +
  theme(axis.text.x = element_text(angle = 0, hjust = .5))
cat("\n=== per-cell core score, by frozen compartment ===\n")
print(s[, .(n = .N, median = round(median(score), 3),
            pct_above_epi_median = round(100 * mean(score > med), 1)), by = froz])
cat("\n=== the same, by cluster ===\n")
print(s[, .(froz = froz[1], median = round(median(score), 3),
            pct_above_epi_median = round(100 * mean(score > med), 1)),
        by = cluster][order(-median)])

fig <- p_a / p_b / p_c + plot_layout(heights = c(1, 2.5, 2.5))
W <- 8.3; H <- 10
for (ext in c("pdf", "png")) {
  f <- file.path(OUT, sprintf("%s_res%s_epi_markers_check.%s", SAMPLE, RES, ext))
  if (ext == "png" && requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(f, width = W, height = H, units = "in", res = 400, background = "white")
  } else if (ext == "png") {
    png(f, width = W, height = H, units = "in", res = 400, bg = "white", type = "cairo")
  } else {
    grDevices::cairo_pdf(f, width = W, height = H)
  }
  print(fig); invisible(dev.off())
  log("wrote ", basename(f))
}
log("done -> ", OUT)
