#!/usr/bin/env Rscript
# One-off check: what are the res0.3 clusters actually made of?
#
# The compartment call only ever asks which of three modules wins (novaseq_04_cluster.R
# scores epithelial / immune / stromal and takes the argmax), so a cluster that is
# neither of the three still gets forced into one. This lists each cluster's top
# markers without reference to that panel, so a lineage the panel does not cover --
# endothelium, myeloid, mast, nerve, adipocyte -- shows up as itself.
#
# Seurat's own FindAllMarkers, Wilcoxon, 2,000 cells per cluster. The cap only affects
# power, and with 641 to 210,750 cells per cluster the uncapped test would be driven
# by cluster size rather than by effect. Nothing here is written back to the pipeline.
#
# TWO STEPS, because neither env can do both halves:
#   - stcnv-r-tuscan has the working Seurat (4.4.0) but SeuratObject 4.1.4 cannot
#     address a v5 Assay5 object at all: ncol() and colnames() come back empty,
#     Idents() falls through to spatial_snn_res.1, and every comparison fails with
#     "Cannot find cells provided";
#   - stcnv-r-shared reads the object fine but pairs Seurat 5.1.0 with SeuratObject
#     5.4.0, which made the `slot` argument of GetAssayData() defunct, so
#     FindAllMarkers stops with a lifecycle error on every comparison and returns
#     nothing.
# Step 1 opens the object in stcnv-r-shared and writes out the capped data matrix
# plus idents; step 2 rebuilds a v4 object from those in stcnv-r-tuscan and runs
# Seurat's FindAllMarkers on it. Capping before the call rather than through
# max.cells.per.ident is equivalent -- that argument downsamples every ident too.
#
# Run (the resolution defaults to 0.3; pass another to repeat the whole thing there):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05d_cluster_markers.R extract [res]
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05d_cluster_markers.R test [res]

suppressWarnings(suppressMessages({ library(Seurat); library(data.table); library(Matrix) }))

ARGS <- commandArgs(trailingOnly = TRUE)
STEP <- ARGS[1]
if (is.na(STEP) || !STEP %in% c("extract", "test"))
  stop("usage: novaseq_05d_cluster_markers.R extract|test [res]")
RES  <- if (length(ARGS) > 1) ARGS[2] else "0.3"

PROJ   <- "<DATA_ROOT>/high_res_ST"
SAMPLE <- "IDC_Grade3"
CLUST  <- file.path(PROJ, "novaseq/processed_data/clustering", SAMPLE)
OUT    <- file.path(PROJ, "results", "strataMap_module_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
BRIDGE <- file.path(OUT, sprintf("%s_res%s_markers_input.rds", SAMPLE, RES))

TOP_N <- 10L
CAP   <- 2000L
PANEL <- c(epithelial = "EPCAM,KRT8,KRT18,KRT19,MUC1",
           immune     = "CD3D,CD3E,JCHAIN,MZB1,MS4A2",
           stromal    = "COL1A1,COL1A2,LUM,THBS2,VWF,PECAM1")
panel_of <- unlist(lapply(names(PANEL), function(k)
  setNames(rep(k, length(strsplit(PANEL[[k]], ",")[[1]])), strsplit(PANEL[[k]], ",")[[1]])))

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## ============================================================ step 1: extract
if (STEP == "extract") {
  stopifnot(packageVersion("SeuratObject") >= "5.0.0")
  # only the res0.3 annotated object exists, but its metadata carries every
  # resolution's spatial_snn column, so it is the source whatever RES is
  log("reading the annotated object (7 GB gz, a few minutes)")
  seu <- readRDS(file.path(CLUST, paste0(SAMPLE, "_seurat_annotated_res0.3.rds.gz")))
  # res 0.3 is the frozen call and carries cluster_frozen; any other resolution is
  # read straight off its spatial_snn column
  ic  <- if (RES == "0.3") "cluster_frozen" else paste0("spatial_snn_res.", RES)
  if (!ic %in% colnames(seu@meta.data)) stop("no such ident column: ", ic)
  cl  <- as.character(seu@meta.data[[ic]]); names(cl) <- colnames(seu)
  log("cells: ", ncol(seu), "  genes: ", nrow(seu), "  clusters: ", length(unique(cl)))

  set.seed(1)
  keep <- unlist(lapply(split(names(cl), cl),
                        function(i) if (length(i) > CAP) sample(i, CAP) else i),
                 use.names = FALSE)
  # LayerData is the v5 accessor; GetAssayData(slot=) is the one that is defunct here
  X <- LayerData(seu, assay = "spatial", layer = "data")[, keep, drop = FALSE]
  saveRDS(list(data = X, idents = cl[keep]), BRIDGE)
  log("wrote ", BRIDGE, "  (", nrow(X), " genes x ", ncol(X), " cells)")
  log("now run step 'test' in stcnv-r-tuscan")
  quit(save = "no")
}

## =============================================================== step 2: test
stopifnot(packageVersion("SeuratObject") < "5.0.0")
b <- readRDS(BRIDGE)
log("rebuilding a v4 object: ", nrow(b$data), " genes x ", ncol(b$data), " cells")
seu <- CreateSeuratObject(counts = b$data, assay = "spatial", min.cells = 0, min.features = 0)
seu <- SetAssayData(seu, assay = "spatial", slot = "data", new.data = b$data)
Idents(seu) <- factor(b$idents[colnames(seu)],
                      levels = as.character(sort(as.integer(unique(b$idents)))))
stopifnot(!anyNA(Idents(seu)), ncol(seu) == length(b$idents))
log("clusters: ", paste(levels(Idents(seu)), collapse = ","))

log("FindAllMarkers (Seurat ", as.character(packageVersion("Seurat")), ", wilcox)")
m <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25,
                    test.use = "wilcox", verbose = FALSE)
M <- as.data.table(m); M[, cluster := as.character(cluster)]
fwrite(M, file.path(OUT, sprintf("%s_res%s_findallmarkers.csv", SAMPLE, RES)))
log("markers passing the thresholds: ", nrow(M))

# top N per cluster by effect size among the significant ones
top <- M[p_val_adj < 0.05][order(cluster, -avg_log2FC), head(.SD, TOP_N), by = cluster]
top[, rank := seq_len(.N), by = cluster]
top[, panel := fifelse(gene %in% names(panel_of), panel_of[gene], "")]

af <- file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES))
if (file.exists(af)) {                       # only the frozen resolution has one
  ann <- fread(af)[, .(cluster = as.character(cluster), n_cells,
                       call = major_celltype, margin, ambiguous)]
} else {
  cr  <- fread(file.path(CLUST, sprintf("%s_compartment_by_resolution.csv", SAMPLE)))
  ann <- cr[as.character(resolution) == RES,
            .(cluster = as.character(cluster), n_cells, call = compartment,
              margin, ambiguous = margin < 0.5)]
}
top[ann, `:=`(n_cells = i.n_cells, call = i.call,
              margin = i.margin, ambiguous = i.ambiguous), on = "cluster"]
setcolorder(top, c("cluster", "n_cells", "call", "ambiguous", "rank", "gene",
                   "panel", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
fwrite(top, file.path(OUT, sprintf("%s_res%s_top%d_markers.csv", SAMPLE, RES, TOP_N)))

# readable summary: one line per cluster, panel members marked
cat("\n")
for (cl in ann[order(-n_cells)]$cluster) {
  t <- top[cluster == cl]
  if (!nrow(t)) { cat(sprintf("%-3s n=%7d  %-6s  (no markers passed)\n", cl,
                              ann$n_cells[ann$cluster == cl],
                              ann$call[ann$cluster == cl])); next }
  tag <- c(epithelial = "^E", immune = "^I", stromal = "^S")
  g <- paste0(t$gene, ifelse(t$panel == "", "", tag[t$panel]))
  cat(sprintf("%-3s n=%7d  %-6s%s  %s\n", cl, t$n_cells[1], t$call[1],
              ifelse(t$ambiguous[1], "*", " "), paste(g, collapse = ", ")))
}
cat("\n^E/^I/^S = already in the epithelial/immune/stromal panel; unmarked genes are not.\n")
log("done -> ", OUT)
