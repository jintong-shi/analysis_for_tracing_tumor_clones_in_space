suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(stringr)
  #library(stringdist)
})
## FFPE 
DEFAULT_DATASET_DIR  <- "<DATA_ROOT>/realdata/datasets/Ding_nature_wes_ff_nobam"
DEFAULT_OUT_BASENAME <- "wes_segments_results.csv"
DEFAULT_SAMPLE_COL   <- "sample"
DATASET_WES_SUBDIR        <- "truth_links/wes_cnv"
DATASET_WES_PATTERN   <- "_wes_cnv.csv"
DATASET_FLAT_RUN_SUBDIR   <- "run_links"
DATASET_DISCARD_PREFIX    <- "discard_sample_"
DATASET_SAMPLE_TABLE      <- "samples/sample_table.csv"
DATASET_OUT_SUBDIR        <- "evaluation"
DATASET_RUN_COL           <- "sample"
DATASET_WES_ID_COL        <- "wes_sample_id"
RAW_COLS <- c("run_sample", "patient_id", "CONTIG", "START", "END",
              "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO", "CALL", "sample")

wes_dir      <- file.path(DEFAULT_DATASET_DIR, DATASET_WES_SUBDIR)
sample_table <- file.path(DEFAULT_DATASET_DIR, DATASET_SAMPLE_TABLE)
out_dir      <- file.path(DEFAULT_DATASET_DIR, DATASET_OUT_SUBDIR)
sample_col <- DEFAULT_SAMPLE_COL

wes_dt <- read_xlsx('<DATA_ROOT>/realdata/datasets/Ding_nature_wes_ff_nobam/truth_links/wes_cnv/41586_2024_8087_MOESM4_ESM.xlsx',sheet = 'WES_CNV')
runlink_to_wes <- list(
  # Ding_nature_wes_ff_nobam
  "HT231P1-S1H3Fc2U1Z1Bs1-H3Bs2-Test" = "HT231P1-S1H3Y1D1_1",
  "HT260C1-Th1K1Fc2U1Z1Bs1"            = "HT260C1-Th1K1A3D1_1-13",
  "HT263B1-S1H1Fc2U1Z1Bs1"             = "HT263B1-S1H1A3D1_1",
  "HT268B1-Th1H3Fc2U12Z1Bs1"           = "HT268B1-Th1H3Y1D1_1",
  "HT268B1-Th1H3Fc2U22Z1Bs1"           = "HT268B1-Th1H3Y1D1_1",
  "HT268B1-Th1H3Fc2U2Z1Bs1"            = "HT268B1-Th1H3Y1D1_1",
  "HT271B1-S1H3Fc2U1Z1Bs1"             = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U2Z1Bs1"             = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U3Z1Bs1"             = "HT271B1-S1H3Y1D1_1",
  "HT271B1-S1H3Fc2U4Z1Bs1"             = "HT271B1-S1H3Y1D1_1",
  # Ding_nature_wes_ffpe
  "HT270P1-S1H1Fs5U1Bp1"               = "HT270P1-S1H1A4K1D1_1",
  "HT397B1-H2U1"                       = "merged_HT397B1-S1H4A4D1_1-12",
  "HT397B1-H3U1"                       = "merged_HT397B1-S1H4A4D1_1-12",
  "HT397B1-H3U2"                       = "merged_HT397B1-S1H4A4D1_1-12",
  #"discard_sample_HT397B1-H3U4"        = "HT397B1-W1D1_1-15",
  #"discard_sample_HT397B1-H3U6"        = "HT397B1-W1D1_1-15",
  "HT425"                              = "merged_HT425B1-S1H1A3D1_1-12",
  "HT486B1-S1H1Fs1U1Bp1"               = "HT486B1-S1H1A3Y1D1_1"
)
for (sample in names(runlink_to_wes)){
  wes_dt_sub <- dplyr::filter(wes_dt,sample==sample)
  write.csv(wes_dt_sub,paste0('<DATA_ROOT>/realdata/datasets/Ding_nature_wes/truth_links/wes_cnv/',sample,"_wes_cnv.csv"))
}

