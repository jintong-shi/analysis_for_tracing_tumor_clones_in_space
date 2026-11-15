#!/usr/bin/env Rscript
# =====================================================================
# ONE-TIME input_data builder (provenance step, NOT a figure script).
#
# Fig 5 panel a is a schematic of ONE representative tissue patch. This script
# performs the (deterministic) window curation: it scans the segmentation grid
# for the most compact well-mixed epithelial/other meta-cell window (identical
# selection logic to tests/schematic_units/plot_units_schematic.R -> focus=164),
# then writes ONLY the raw 2um bins inside that window (x, y, object_id, is_epi,
# meta pseudo_id) plus the window bounds. All schematic geometry (8/16um bins,
# cell outlines, meta-cell hulls, palette) is recomputed later in
# scripts/fig5_main.R from these raw bins.
#
# Run under stcnv-r-tuscan:
#   <CONDA_ROOT>/envs/stcnv-r-tuscan/bin/Rscript \
#     paper/figures_final/high_res_st/input_data/prep/prep_schematic_window.R
# =====================================================================
suppressWarnings(suppressMessages(library(data.table)))
base   <- "<DATA_ROOT>/high_res_ST/visiumHD/processed_data"
map_f  <- file.path(base, "deconvolution/processed_data/histo_decon_pure_epithelial_cell_mapping_coord.csv")
seg_f  <- file.path(base, "segmentation/metadata_labels_joint.csv")
meta_f <- file.path(base, "pseudo_spot/meta_r40_Smax20.csv")
OD     <- "../schematic"
dir.create(OD, showWarnings = FALSE, recursive = TRUE)

## epithelial reference + meta-cell membership ------------------------------
map  <- fread(map_f)
meta <- fread(meta_f); setnames(meta, c("object_id", "pseudo_id"))
epi_set <- unique(map$object_id)
map <- merge(map, meta, by = "object_id", all.x = TRUE)

## load the whole tissue-cell mask cropped to the epi region (awk, no full read) --
crop <- 60L
bx <- c(min(map$x) - crop, max(map$x) + crop); by <- c(min(map$y) - crop, max(map$y) + crop)
awkcmd <- sprintf("awk -F, 'NR>1 && $4>=%d && $4<=%d && $3>=%d && $3<=%d && $11!=0 {print $4\",\"$3\",\"$11}' %s",
                  bx[1], bx[2], by[1], by[2], seg_f)
tis <- fread(cmd = awkcmd, header = FALSE, col.names = c("x", "y", "oid"))
tis[, is_epi := oid %in% epi_set]

## pick the compact well-mixed window (identical to plot_units_schematic.R) -----
half <- 16L; hb <- 5L; nb <- (2 * hb + 1)^2
b8 <- tis[, .(n = .N, fe = mean(is_epi)), by = .(x8 = x %/% 4, y8 = y %/% 4)][n >= 8]
b8[, type := ifelse(fe >= 0.5, 1L, -1L)]
ox <- min(b8$x8); oy <- min(b8$y8)
gx <- max(b8$x8) - ox + 1; gy <- max(b8$y8) - oy + 1
TY <- matrix(0L, gx, gy); TY[cbind(b8$x8 - ox + 1, b8$y8 - oy + 1)] <- b8$type
pad <- matrix(0L, gx + 2, gy + 2); pad[2:(gx + 1), 2:(gy + 1)] <- TY
core <- pad[2:(gx + 1), 2:(gy + 1)]
U <- (core * pad[2:(gx + 1), 3:(gy + 2)] < 0) + (core * pad[2:(gx + 1), 1:gy] < 0) +
     (core * pad[1:gx, 2:(gy + 1)] < 0) + (core * pad[3:(gx + 2), 2:(gy + 1)] < 0)
intimg <- function(M) { IP <- t(apply(apply(M, 2, cumsum), 1, cumsum)); rbind(0, cbind(0, IP)) }
IU <- intimg(U); IE <- intimg(TY == 1); IO <- intimg(TY == -1)
wsum <- function(IP, mx, my) {
  x0 <- max(mx - hb, 1); x1 <- min(mx + hb, gx); y0 <- max(my - hb, 1); y1 <- min(my + hb, gy)
  IP[x1 + 1, y1 + 1] - IP[x0, y1 + 1] - IP[x1 + 1, y0] + IP[x0, y0] }
cand <- map[!is.na(pseudo_id), .(cx = round(mean(x)), cy = round(mean(y)),
                                 ncell = uniqueN(object_id),
                                 spread = sqrt(mean((x - mean(x))^2 + (y - mean(y))^2))),
             by = pseudo_id][ncell >= 10 & ncell <= 14]
cand[, `:=`(mx = cx %/% 4 - ox + 1, my = cy %/% 4 - oy + 1)]
cand <- cand[mx - hb >= 1 & mx + hb <= gx & my - hb >= 1 & my + hb <= gy]
cand[, `:=`(mix = mapply(function(a, b) wsum(IU, a, b), mx, my),
            ne  = mapply(function(a, b) wsum(IE, a, b), mx, my),
            no  = mapply(function(a, b) wsum(IO, a, b), mx, my))]
cand <- cand[ne >= 25 & no >= 25 & (ne + no) >= 0.75 * nb]
focus <- cand[which.min(spread)]
xr <- c(focus$cx - half, focus$cx + half); yr <- c(focus$cy - half, focus$cy + half)
focus_id <- focus$pseudo_id

## every in-tissue 2um bin in the chosen window (incl. unassigned oid=0) ---------
awk2 <- sprintf("awk -F, 'NR>1 && $4>=%d && $4<=%d && $3>=%d && $3<=%d {print $4\",\"$3\",\"$11}' %s",
                xr[1], xr[2], yr[1], yr[2], seg_f)
wall <- fread(cmd = awk2, header = FALSE, col.names = c("x", "y", "oid"))
wall[is.na(oid), oid := 0L]
wall[, is_epi := oid %in% epi_set]
wall <- merge(wall, meta, by.x = "oid", by.y = "object_id", all.x = TRUE)   # pseudo_id (NA if none)
setcolorder(wall, c("x", "y", "oid", "is_epi", "pseudo_id"))
setorder(wall, x, y)

fwrite(wall, file.path(OD, "schematic_window_bins.csv"))
fwrite(data.table(xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2], focus = focus_id),
       file.path(OD, "schematic_window_bounds.csv"))
cat(sprintf("focus=%s window %dx%d bins=%d (labelled=%d) epi_cells=%d\n",
            focus_id, diff(xr) + 1, diff(yr) + 1, nrow(wall), sum(wall$oid != 0),
            uniqueN(wall[is_epi == TRUE & oid != 0]$oid)))
