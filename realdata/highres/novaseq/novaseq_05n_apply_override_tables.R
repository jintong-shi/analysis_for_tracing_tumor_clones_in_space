#!/usr/bin/env Rscript
# Apply the novaseq_05_annotate.R manual override (IDC_Grade3 res0.3 cluster 7
# imm -> stroma) to the two annotation CSVs only, so figures can use the new call without
# the ~25 min annotation rerun (most of it is rewriting the 7 GB annotated rds).
# Writes exactly the major_celltype / override columns 05_annotate would; the rds and the
# _annot_res0.3_* plots stay stale until that rerun. Originals are kept as *.pre_override.
# Idempotent: a second run changes nothing.
#
# Rationale: docs/IDC_Grade3_compartment_override.md
#
# Run (~30 s):
#   <CONDA_ROOT>/envs/stcnv-r-shared/bin/Rscript \
#     high_res_ST/novaseq/code/novaseq_05n_apply_override_tables.R

suppressWarnings(suppressMessages(library(data.table)))

SAMPLE <- "IDC_Grade3"; RES <- "0.3"
OVERRIDE <- c(`7` = "stroma")          # same as novaseq_05_annotate.R
CLUST <- file.path("<DATA_ROOT>/high_res_ST",
                   "novaseq", "processed_data", "clustering", SAMPLE)
f_tab   <- file.path(CLUST, sprintf("%s_annotation_res%s.csv", SAMPLE, RES))
f_cells <- file.path(CLUST, sprintf("%s_annotated_cells_res%s.csv", SAMPLE, RES))
for (f in c(f_tab, f_cells)) {
  b <- paste0(f, ".pre_override")
  if (!file.exists(b)) stopifnot(file.copy(f, b))
}

tab <- fread(f_tab, colClasses = list(character = "cluster"))
epi_before <- tab[major_celltype == "epi", sum(n_cells)]
if (!"override" %in% names(tab)) tab[, override := ""]
i <- match(names(OVERRIDE), tab$cluster)
tab[i, override := ifelse(major_celltype == OVERRIDE, override, paste0(major_celltype, "->", OVERRIDE))]
tab[i, major_celltype := unname(OVERRIDE)]
setcolorder(tab, c("cluster", "n_cells", "major_celltype", "compartment", "margin",
                   "ambiguous", "override"))
stopifnot(tab[major_celltype == "epi", sum(n_cells)] == epi_before)
fwrite(tab, f_tab)
print(tab[, .(cluster, n_cells, major_celltype, override)])

cells <- fread(f_cells)
if (names(cells)[1] == "V1") setnames(cells, 1, "")   # keep write.csv's unnamed row-name column
hit <- as.character(cells$cluster_frozen) %in% names(OVERRIDE)
cells[hit, major_celltype := OVERRIDE[as.character(cluster_frozen)]]
fwrite(cells, f_cells)
print(cells[, .N, by = major_celltype])
