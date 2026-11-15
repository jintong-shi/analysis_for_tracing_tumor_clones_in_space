#!/usr/bin/env Rscript
# Per-cell immunoglobulin counts for all 715,479 cells -- the raw material for the
# ambient-Ig figure.
#
# The question is whether the Ig-high cells are plasma cells or cells sitting in an
# Ig-rich neighbourhood. That cannot be answered from a normalised matrix capped at
# 2,000 cells per cluster, so this reads the filtered h5 directly.
#
# One pass over matrix/indices + matrix/data (1.23 billion non-zeros) accumulates every
# Ig group at once -- the cost is reading the two arrays, not the number of genes, so
# counting 186 genes costs what counting one gene costs. Columns are recovered from
# indptr with findInterval, exactly as prep_strataMap_extracts.R does for EPCAM.
#
# Groups are kept apart because they answer different questions:
#   IGH_C  isotype constant regions   -- how much Ig protein message is in the cell
#   IGK    kappa (IGKC + IGKV)        -- light-chain identity, see below
#   IGL    lambda (IGLC + IGLV)       --   "
#   IGHV   heavy variable regions     -- clonality
#   JCHAIN, MZB1                      -- plasma-cell machinery, NOT secreted product
#
# A real plasma cell is light-chain restricted: it transcribes kappa or lambda, never
# both. Ambient mRNA is a pool drawn from every plasma cell in the neighbourhood, so a
# contaminated non-plasma cell carries the local kappa:lambda MIXTURE. That is the one
# test ambient cannot fake, and it needs the two light chains counted separately.
# MZB1/JCHAIN are included as the control: they are made by plasma cells but not
# secreted in bulk, so they contaminate far less than IGKC does.
#
# Run (~10 min, mostly the h5 scan):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05i_ig_ambient_extract.R

suppressWarnings(suppressMessages({ library(hdf5r); library(data.table) }))

PROJ   <- "<DATA_ROOT>"
SAMPLE <- "IDC_Grade3"
H5     <- file.path(PROJ, "high_res_ST/novaseq/processed_data/cells",
                    paste0(SAMPLE, "_cells_filtered.h5"))
OUT    <- file.path(PROJ, "high_res_ST/results", "strataMap_module_scores")
CHUNK  <- 5e7

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

h5 <- H5File$new(H5, "r")
g  <- h5[["matrix/features/name"]]$read()
bc <- h5[["matrix/barcodes"]]$read()

# explicit constant-region list: ^IGH[AGMDE] also matches IGHMBP2 and IGHEP2, which are
# not immunoglobulins at all
IGH_C <- c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHM","IGHD","IGHE")
GRP <- list(IGH_C = intersect(IGH_C, g),
            IGHV  = grep("^IGHV", g, value = TRUE),
            IGK   = grep("^IGK[CV]", g, value = TRUE),
            IGL   = grep("^IGL[CV]", g, value = TRUE),
            JCHAIN = "JCHAIN", MZB1 = "MZB1")
stopifnot(all(lengths(GRP) > 0))
log(paste(sprintf("%s=%d", names(GRP), lengths(GRP)), collapse = "  "))

# row (0-based) -> group index; 0 means "not an Ig gene"
grp_of <- integer(length(g))
for (k in seq_along(GRP)) grp_of[match(GRP[[k]], g)] <- k

ptr <- h5[["matrix/indptr"]]$read()                    # 0-based, length ncol+1
di  <- h5[["matrix/indices"]]; dd <- h5[["matrix/data"]]
nnz <- di$dims
ncell <- length(bc)
M <- matrix(0, nrow = ncell, ncol = length(GRP), dimnames = list(NULL, names(GRP)))
total <- numeric(ncell)                                 # every cell's UMI total

log("scanning ", format(nnz, big.mark = ","), " non-zeros in ",
    ceiling(nnz / CHUNK), " chunks")
for (st in seq(1, nnz, by = CHUNK)) {
  en  <- min(st + CHUNK - 1, nnz)
  idx <- di$read(args = list(st:en))                    # 0-based gene rows
  val <- as.numeric(dd$read(args = list(st:en)))
  col <- findInterval(seq.int(st, en) - 1L, ptr)        # 1-based cell
  s <- rowsum(val, col); r <- as.integer(rownames(s))
  total[r] <- total[r] + s[, 1]
  gk <- grp_of[idx + 1L]
  for (k in seq_along(GRP)) {
    sel <- gk == k
    if (!any(sel)) next
    s <- rowsum(val[sel], col[sel]); r <- as.integer(rownames(s))
    M[r, k] <- M[r, k] + s[, 1]
  }
  log(sprintf("  %3.0f%%", 100 * en / nnz))
}
h5$close_all()

d <- data.table(cell = bc, nCount = total)
for (k in names(GRP)) d[[k]] <- M[, k]
f <- file.path(OUT, sprintf("%s_percell_ig_counts.csv", SAMPLE))
fwrite(d, f)
log("wrote ", f, "  (", nrow(d), " cells)")
log(sprintf("median Ig fraction (IGH_C+IGK+IGL): %.3f%%",
            100 * median((d$IGH_C + d$IGK + d$IGL) / d$nCount)))
