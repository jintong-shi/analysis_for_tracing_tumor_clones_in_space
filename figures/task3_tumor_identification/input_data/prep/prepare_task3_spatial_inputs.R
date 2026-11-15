#!/usr/bin/env Rscript
# =====================================================================
# FREEZE STEP (run once): build self-contained spatial inputs for Fig 4.
# Reads pipeline outputs (predictions_standardized + tissue_positions) and
# writes one long CSV per sample (spatial_<stem>.csv) into ../real. The figure
# scripts read ONLY these frozen CSVs (they never reach back into the pipeline):
# fig4_task3_main.R (V1_2/H2_2/H1_5) and figS_task3_supp.R (P4_P13/pt15_P2).
#
# Per-spot class (unified ids; see figure_settings.tsv spatial_truth/spatial_pred):
#   truth  panel : truth_bin 1->pos, 0->neg, else->other
#   pred   panel : pred_bin  1->pos, 0->neg, present-but-blank->unassigned,
#                            spot not in tool file->na
#   Numbat never shows na/unassigned -> those become neg (diploid). On samples
#   where Numbat was not run (V1_2) the whole Numbat panel is diploid.
#
# 5 tools = CopyKAT/SCEVAN/CopyVAE/TUSCAN @ AllSpots_no-ref + Numbat @
# Decon_StromaImm-ref (FF); FFPE samples have no Numbat -> 4 tools.
# =====================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

# Run this from the prep/ directory. figure_setup.R lives in figures/utils/;
# outputs are written to ../real (figures/task3_tumor_identification/input_data/real),
# which fig4_task3_main.R and figS_task3_supp.R read as IN_REAL.
source("../../../utils/figure_setup.R")
PRED <- "<DATA_ROOT>/realdata_unified_by_dataset/results/task3_tumor_vs_histology/predictions_standardized"
DATS <- "<DATA_ROOT>/realdata_unified_by_dataset/datasets"
OUT  <- "../real"

read_positions <- function(path) {
  cols <- c("barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres")
  has_header <- grepl("barcode", readLines(path, n = 1), ignore.case = TRUE)
  d <- if (has_header) read_csv(path, show_col_types = FALSE)
       else            read_csv(path, col_names = cols, show_col_types = FALSE)
  d %>% filter(in_tissue == 1) %>%
    transmute(spot_id = barcode, x = pxl_col_in_fullres, y = pxl_row_in_fullres)
}
std_csv <- function(ds, s, run, tool)
  file.path(PRED, ds, s, run, tool, "tumor_normal_standardized.csv")

build_long <- function(ds, s, pos, panels, truth_from) {
  universe <- read_positions(pos)
  panel_levels <- c("Histology", vapply(panels, `[[`, "", "label"))

  truth <- read_csv(std_csv(ds, s, truth_from$run, truth_from$tool),
                    show_col_types = FALSE) %>%
    transmute(spot_id, panel = "Histology", panel_role = "truth",
              class = case_when(truth_bin == 1 ~ "pos",
                                truth_bin == 0 ~ "neg", TRUE ~ "other"))
  preds <- lapply(panels, function(p) {
    csv <- std_csv(ds, s, p$run, p$tool)
    if (!file.exists(csv))
      return(tibble(spot_id = character(), panel = p$label,
                    panel_role = "pred", class = character()))
    read_csv(csv, show_col_types = FALSE) %>%
      transmute(spot_id, panel = p$label, panel_role = "pred",
                class = case_when(pred_bin == 1 ~ "pos",
                                  pred_bin == 0 ~ "neg", TRUE ~ "unassigned"))
  })
  assigned <- bind_rows(truth, bind_rows(preds))
  numbat_lab <- vapply(panels, function(p)
    if (identical(p$tool, "numbat")) p$label else NA_character_, "")
  numbat_lab <- numbat_lab[!is.na(numbat_lab)]

  roles <- c(Histology = "truth",
             setNames(rep("pred", length(panels)),
                      vapply(panels, `[[`, "", "label")))
  crossing(universe, panel = panel_levels) %>%
    left_join(assigned %>% select(spot_id, panel, class),
              by = c("spot_id", "panel")) %>%
    mutate(panel_role = unname(roles[panel]),
           class = ifelse(is.na(class), ifelse(panel_role == "truth", "other", "na"), class),
           class = ifelse(panel %in% numbat_lab & class %in% c("na", "unassigned"),
                          "neg", class),
           panel_order = match(panel, panel_levels),
           sample = s) %>%
    select(sample, panel, panel_order, panel_role, spot_id, x, y, class) %>%
    arrange(panel_order)
}

mk <- function(tool, run) list(run = run, tool = tool, label = unname(TOOL_DISPLAY[tool]))
EXPR <- lapply(c("copykat", "scevan", "copyvae", "tuscan"),
               function(t) mk(t, "AllSpots_no-ref"))
FIVE <- c(EXPR, list(mk("numbat", "Decon_StromaImm-ref")))
FOUR <- EXPR
TF   <- list(run = "AllSpots_no-ref", tool = "copyvae")
pstar <- function(ds, root, s)
  file.path(DATS, ds, root, s, "AllSpots_no-ref",
            "STARCH/starch_10x_subset/tissue_positions_list.csv")

jobs <- list(
  list(stem="V1_2",    ds="nature_prostate_ff",       pos=pstar("nature_prostate_ff","run_links","V1_2"),  pn=FIVE),
  list(stem="H2_2",    ds="nature_prostate_ff",       pos=pstar("nature_prostate_ff","run_links","H2_2"),  pn=FIVE),
  list(stem="H1_5",    ds="nature_prostate_ff",       pos=pstar("nature_prostate_ff","run_links","H1_5"),  pn=FIVE),
  list(stem="P4_P13",  ds="prostate_eva_ffpe",        pos=pstar("prostate_eva_ffpe","run_CNV","P4_P13"),   pn=FOUR),
  list(stem="pt15_P2", ds="metastatic_prostate_ffpe",
       pos=file.path(DATS,"metastatic_prostate_ffpe/source_links",
                     "selected_oxford_orig_data/pt15_P2/spatial/tissue_positions.csv"), pn=FOUR))

for (j in jobs) {
  d <- build_long(j$ds, j$stem, j$pos, j$pn, TF)
  f <- file.path(OUT, paste0("spatial_", j$stem, ".csv"))
  write_csv(d, f)
  cat(sprintf("[freeze] %s  (%d spots x %d panels)\n",
              basename(f), length(unique(d$spot_id)), length(unique(d$panel))))
}
