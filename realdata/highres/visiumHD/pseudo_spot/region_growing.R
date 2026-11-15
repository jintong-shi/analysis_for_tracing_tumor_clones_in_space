## VisiumHD AP5 — epithelial meta-cells by region growing (step 1/2)
##
## Groups the bin2cell/RCTD pure-epithelial single cells into spatially compact
## "meta-cells" (pseudo-spots): connect every pair of cells within r_adj µm,
## take connected components, then BFS-split any component larger than S_max so
## no meta-cell exceeds that size. Components smaller than S_min are dropped
## (their cells keep pseudo_id = NA and are lost downstream).
##
## Outputs, per parameter pair:
##   meta_r<r>_Smax<S>.csv          cell (object_id) -> pseudo_id, NA if dropped
##   epi_pseudo_r<r>_Smax<S>.rds    the same Seurat object with $pseudo_id
## Step 2/2 (generate_pseudo_spot_seu.R) aggregates these into the pseudo-spot
## object InferCNV is run on.
##
## Re-running this OVERWRITES meta_r40_Smax20.csv. The copy currently on disk is
## from an older run (65,278 cells -> 5,614 meta-cells) and does NOT correspond to
## epi_pseudo_r40_Smax20*.rds (54,098 -> 4,721): the two share pseudo_id numbering
## but not its meaning and must never be joined (see the note in
## figures/high_res_st/input_data/prep/prep_subclone_coords.R). Re-running makes
## the CSV and the object consistent again.
##
## Run under stcnv-r-shared:
##   Rscript region_growing.R
suppressMessages({
  library(Seurat)
  library(dbscan)   # frNN
  library(igraph)
})

## paths -------------------------------------------------------------------
epi_rds <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/deconvolution/processed_data/histo_decon_pure_epi_seu.rds.gz"
outdir  <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data/pseudo_spot"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## parameters ---------------------------------------------------------------
## The paper uses r_adj = 40 µm, S_max = 20. The original scan was
## r_adj 25/30/35/40 x S_max 20/30/40/50 — widen the vectors to reproduce it.
r_adj_vals <- c(40)
S_max_vals <- c(20)
S_min      <- 3

epi.seu.hd <- readRDS(epi_rds)
coords     <- epi.seu.hd@reductions$spatial@cell.embeddings

## region growing -----------------------------------------------------------
region_grow <- function(coords, r_adj, S_max) {
  nbr  <- frNN(coords, eps = r_adj, sort = FALSE)
  lens <- lengths(nbr$id)
  edges_df <- data.frame(from = rep(seq_along(nbr$id), times = lens),
                         to   = unlist(nbr$id))
  edges_df <- unique(edges_df[edges_df$from != edges_df$to, ])
  g <- graph_from_data_frame(edges_df, directed = FALSE,
                             vertices = seq_len(nrow(coords)))
  comp <- igraph::components(g)$membership

  split_by_size <- function(idx) {
    if (length(idx) <= S_max) return(list(idx))
    subg   <- induced_subgraph(g, idx)
    unseen <- V(subg)
    subs   <- list()
    while (length(unseen) > 0) {
      hop <- subcomponent(subg, unseen[1], mode = "all")
      if (length(hop) > S_max) hop <- hop[1:S_max]
      subs[[length(subs) + 1]] <- idx[hop]
      unseen <- setdiff(unseen, hop)
    }
    subs
  }
  raw <- unlist(lapply(split(seq_along(comp), comp), split_by_size),
                recursive = FALSE)
  lapply(raw, unique)
}

for (r_adj in r_adj_vals) {
  for (S_max in S_max_vals) {
    cat(sprintf("### r_adj=%d, S_max=%d ###\n", r_adj, S_max))
    pseudo_list <- Filter(function(x) length(x) >= S_min,
                          region_grow(coords, r_adj, S_max))

    epi <- epi.seu.hd
    epi$pseudo_id <- NA
    for (i in seq_along(pseudo_list)) epi$pseudo_id[pseudo_list[[i]]] <- i

    stem <- sprintf("r%d_Smax%d", r_adj, S_max)
    write.csv(data.frame(cell = colnames(epi), pseudo_id = epi$pseudo_id),
              file.path(outdir, sprintf("meta_%s.csv", stem)), row.names = FALSE)
    saveRDS(epi, file.path(outdir, sprintf("epi_pseudo_%s.rds", stem)))

    cat(sprintf("meta-cells: %d | cells kept: %d / %d\n", length(pseudo_list),
                sum(!is.na(epi$pseudo_id)), ncol(epi)))
  }
}
