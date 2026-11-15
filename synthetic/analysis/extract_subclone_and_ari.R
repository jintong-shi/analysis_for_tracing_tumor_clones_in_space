#!/usr/bin/env Rscript
# Extract each tool's predicted spot-level subclones and score them against the
# ground truth with ARI, per Rep_* batch.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) stop("Cannot resolve script path.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "subclone_helpers.R"))

read_gt_hard_table <- function(gt_file) {
  if (!file.exists(gt_file)) return(NULL)

  df <- tryCatch(
    fread(gt_file, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!all(c("spot_id", "subclone") %in% colnames(df))) return(NULL)

  out <- df[, c("spot_id", "subclone"), drop = FALSE]
  out$spot_id <- as.character(out$spot_id)
  out$subclone <- as.character(out$subclone)
  out <- out[!is.na(out$spot_id) & nzchar(out$spot_id) & !is.na(out$subclone) & nzchar(out$subclone), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out
}

# ------------------------------------------------------------------------------
# Main (prediction + ARI)
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--root_dir", type = "character", default = "<DATA_ROOT>/synthetic/processed_data/random_batch_BAF/generate_code/Complete_Batch_v2", help = "Root benchmark folder containing Rep_*"),
  make_option("--tool_out_name", type = "character", default = "Tool_Outputs",
              help = "Tool output folder name under each rep [default %default]"),
  make_option("--std_in_name", type = "character", default = "Standardized_Inputs",
              help = "Standardized input folder name under each rep [default %default]"),
  make_option("--ref_suffix", type = "character", default = "",
              help = "Suffix appended to ref-tool dir names SCEVAN/Numbat (e.g. _TumorFreeNormal) [default none]"),
  make_option("--infercnv_dir_name", type = "character", default = "",
              help = "Exact InferCNV output dir name under Tool_Outputs (e.g. InferCNV_noHMM_TumorFreeNormal); empty = auto-resolve [default none]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$root_dir)) {
  stop("Usage: Rscript extract_subclone_and_ari.R --root_dir <PATH>")
}
if (!dir.exists(opt$root_dir)) stop("root_dir does not exist: ", opt$root_dir)

run_tool_out_name <- opt$tool_out_name
run_std_in_name <- opt$std_in_name

rep_dirs <- list.dirs(opt$root_dir, full.names = TRUE, recursive = FALSE)
rep_dirs <- rep_dirs[grepl("^Rep_", basename(rep_dirs))]
if (length(rep_dirs) == 0) stop("No Rep_* folders found under: ", opt$root_dir)

message(sprintf("[Init] Found %d replicate folders.", length(rep_dirs)))

pred_summary_rows <- list()
ari_rows <- list()

for (rep_dir in rep_dirs) {
  rep_name <- basename(rep_dir)
  message(sprintf("\n=== %s ===", rep_name))

  syn_dir <- file.path(rep_dir, CFG$synthetic_name)
  gt_hard_out <- file.path(syn_dir, CFG$spatial_dir, CFG$gt_hard_out_file)
  gt_hard_df <- read_gt_hard_table(gt_hard_out)
  gt_hard_ok <- !is.null(gt_hard_df)

  if (!gt_hard_ok) {
    message("  [GT] Missing or invalid hard ground-truth subclone table; ARI will be marked as missing_ground_truth.")
  }

  tool_root <- file.path(rep_dir, run_tool_out_name)
  if (!dir.exists(tool_root)) {
    message(sprintf("  [Pred] Missing %s", tool_root))
    next
  }

  pred_tumor_spots <- read_tumor_spots_metadata(
    rep_dir = rep_dir,
    std_in_name = run_std_in_name,
    tumor_only = CFG$pred_tumor_only,
    tumor_prop_cutoff = CFG$pred_tumor_prop_cutoff
  )
  if (CFG$pred_tumor_only) {
    if (is.null(pred_tumor_spots)) {
      message("  [Pred] Tumor filter requested but metadata not found; using all spots.")
    } else {
      message(sprintf("  [Pred] Tumor spots loaded: %d (Tumor_prop > %.2f)", length(pred_tumor_spots), CFG$pred_tumor_prop_cutoff))
    }
  }

  # Spot universe used to inject a single-cluster prediction when a tool ran
  # but found no subclonal structure (see monoclone classifiers in helpers).
  mono_spots <- if (!is.null(pred_tumor_spots)) {
    pred_tumor_spots
  } else if (gt_hard_ok) {
    gt_hard_df$spot_id
  } else {
    NULL
  }

  # 1) InferCNV
  if (nzchar(opt$infercnv_dir_name)) {
    cand <- file.path(tool_root, opt$infercnv_dir_name)
    infer_dir <- if (dir.exists(cand)) cand else NULL
  } else {
    infer_dir <- resolve_infercnv_dir(tool_root)
  }
  infer_df <- extract_infercnv_table(infer_dir)
  infer_df <- subset_spots(infer_df, pred_tumor_spots)
  infer_out <- if (!is.null(infer_dir)) file.path(infer_dir, CFG$pred_out_file_name) else NA_character_
  infer_ok <- !is.null(infer_dir) && write_subclone_table(infer_df, infer_out)
  message(sprintf(
    "  [InferCNV] %s",
    ifelse(infer_ok, sprintf("Saved (spots=%d, subclones=%d)", nrow(infer_df), count_subclones(infer_df)), "No usable output")
  ))
  pred_summary_rows[[length(pred_summary_rows) + 1]] <- data.frame(
    rep = rep_name, tool = "InferCNV", input_dir = ifelse(is.null(infer_dir), NA, infer_dir),
    out_file = ifelse(infer_ok, infer_out, NA),
    n_spots = ifelse(is.null(infer_df), 0L, nrow(infer_df)),
    n_subclones = count_subclones(infer_df),
    status = ifelse(infer_ok, "ok", "missing_or_empty"),
    stringsAsFactors = FALSE
  )
  ari_rows[[length(ari_rows) + 1]] <- compute_ari_row(
    rep_name = rep_name, tool = "InferCNV",
    gt_df = gt_hard_df, pred_df = infer_df,
    gt_file = ifelse(gt_hard_ok, gt_hard_out, NA),
    pred_file = ifelse(infer_ok, infer_out, NA),
    min_overlap = CFG$ari_min_overlap
  )

  # 2) SCEVAN
  scevan_dir <- file.path(tool_root, paste0("SCEVAN", opt$ref_suffix))
  scevan_df <- extract_scevan_table(scevan_dir)
  scevan_df <- subset_spots(scevan_df, pred_tumor_spots)
  scevan_kind <- "native"
  if (is.null(scevan_df) || nrow(scevan_df) == 0) {
    scevan_cause <- classify_scevan_nooutput(scevan_dir)
    if (scevan_cause == "monoclone") {
      scevan_df <- subset_spots(make_monoclone_df(mono_spots), pred_tumor_spots)
      scevan_kind <- "injected_monoclone"
    } else if (scevan_cause == "no_tumor") {
      scevan_kind <- "failure_no_tumor"
    } else {
      scevan_kind <- "missing"
    }
  }
  scevan_out <- file.path(scevan_dir, CFG$pred_out_file_name)
  scevan_ok <- write_subclone_table(scevan_df, scevan_out)
  message(sprintf(
    "  [SCEVAN ] %s%s",
    ifelse(scevan_ok, sprintf("Saved (spots=%d, subclones=%d)", nrow(scevan_df), count_subclones(scevan_df)), "No usable output"),
    ifelse(scevan_kind == "native", "", sprintf(" [%s]", scevan_kind))
  ))
  pred_summary_rows[[length(pred_summary_rows) + 1]] <- data.frame(
    rep = rep_name, tool = "SCEVAN", input_dir = ifelse(dir.exists(scevan_dir), scevan_dir, NA),
    out_file = ifelse(scevan_ok, scevan_out, NA),
    n_spots = ifelse(is.null(scevan_df), 0L, nrow(scevan_df)),
    n_subclones = count_subclones(scevan_df),
    pred_kind = scevan_kind,
    status = ifelse(scevan_ok, "ok", "missing_or_empty"),
    stringsAsFactors = FALSE
  )
  scevan_ari <- compute_ari_row(
    rep_name = rep_name, tool = "SCEVAN",
    gt_df = gt_hard_df, pred_df = scevan_df,
    gt_file = ifelse(gt_hard_ok, gt_hard_out, NA),
    pred_file = ifelse(scevan_ok, scevan_out, NA),
    min_overlap = CFG$ari_min_overlap
  )
  scevan_ari$pred_kind <- scevan_kind
  ari_rows[[length(ari_rows) + 1]] <- scevan_ari

  # 3) STARCH
  starch_dir <- file.path(tool_root, "STARCH")
  starch_df <- extract_starch_table(starch_dir, rep_dir)
  starch_df <- subset_spots(starch_df, pred_tumor_spots)
  starch_out <- file.path(starch_dir, CFG$pred_out_file_name)
  starch_ok <- write_subclone_table(starch_df, starch_out)
  message(sprintf(
    "  [STARCH ] %s",
    ifelse(starch_ok, sprintf("Saved (spots=%d, subclones=%d)", nrow(starch_df), count_subclones(starch_df)), "No usable output")
  ))
  pred_summary_rows[[length(pred_summary_rows) + 1]] <- data.frame(
    rep = rep_name, tool = "STARCH", input_dir = ifelse(dir.exists(starch_dir), starch_dir, NA),
    out_file = ifelse(starch_ok, starch_out, NA),
    n_spots = ifelse(is.null(starch_df), 0L, nrow(starch_df)),
    n_subclones = count_subclones(starch_df),
    status = ifelse(starch_ok, "ok", "missing_or_empty"),
    stringsAsFactors = FALSE
  )
  ari_rows[[length(ari_rows) + 1]] <- compute_ari_row(
    rep_name = rep_name, tool = "STARCH",
    gt_df = gt_hard_df, pred_df = starch_df,
    gt_file = ifelse(gt_hard_ok, gt_hard_out, NA),
    pred_file = ifelse(starch_ok, starch_out, NA),
    min_overlap = CFG$ari_min_overlap
  )

  # 4) Numbat
  numbat_dir <- file.path(tool_root, paste0("Numbat", opt$ref_suffix))
  numbat_df <- extract_numbat_table(numbat_dir)
  numbat_df <- subset_spots(numbat_df, pred_tumor_spots)
  numbat_kind <- "native"
  if (is.null(numbat_df) || nrow(numbat_df) == 0) {
    numbat_kind <- classify_numbat_nooutput(numbat_dir)
    numbat_df <- NULL
  }
  numbat_out <- file.path(numbat_dir, CFG$pred_out_file_name)
  numbat_ok <- write_subclone_table(numbat_df, numbat_out)
  message(sprintf(
    "  [Numbat ] %s%s",
    ifelse(numbat_ok, sprintf("Saved (spots=%d, subclones=%d)", nrow(numbat_df), count_subclones(numbat_df)), "No usable output"),
    ifelse(numbat_kind == "native", "", sprintf(" [%s]", numbat_kind))
  ))
  pred_summary_rows[[length(pred_summary_rows) + 1]] <- data.frame(
    rep = rep_name, tool = "Numbat", input_dir = ifelse(dir.exists(numbat_dir), numbat_dir, NA),
    out_file = ifelse(numbat_ok, numbat_out, NA),
    n_spots = ifelse(is.null(numbat_df), 0L, nrow(numbat_df)),
    n_subclones = count_subclones(numbat_df),
    pred_kind = numbat_kind,
    status = ifelse(numbat_ok, "ok", "missing_or_empty"),
    stringsAsFactors = FALSE
  )
  numbat_ari <- compute_ari_row(
    rep_name = rep_name, tool = "Numbat",
    gt_df = gt_hard_df, pred_df = numbat_df,
    gt_file = ifelse(gt_hard_ok, gt_hard_out, NA),
    pred_file = ifelse(numbat_ok, numbat_out, NA),
    min_overlap = CFG$ari_min_overlap
  )
  numbat_ari$pred_kind <- numbat_kind
  ari_rows[[length(ari_rows) + 1]] <- numbat_ari

  # 5) CalicoST
  calicost_dir <- file.path(tool_root, "CalicoST")
  calicost_df <- extract_calicost_table(calicost_dir)
  calicost_df <- subset_spots(calicost_df, pred_tumor_spots)
  calicost_out <- file.path(calicost_dir, CFG$pred_out_file_name)
  calicost_ok <- write_subclone_table(calicost_df, calicost_out)
  message(sprintf(
    "  [CalicoST] %s",
    ifelse(calicost_ok, sprintf("Saved (spots=%d, subclones=%d)", nrow(calicost_df), count_subclones(calicost_df)), "No usable output")
  ))
  pred_summary_rows[[length(pred_summary_rows) + 1]] <- data.frame(
    rep = rep_name, tool = "CalicoST", input_dir = ifelse(dir.exists(calicost_dir), calicost_dir, NA),
    out_file = ifelse(calicost_ok, calicost_out, NA),
    n_spots = ifelse(is.null(calicost_df), 0L, nrow(calicost_df)),
    n_subclones = count_subclones(calicost_df),
    status = ifelse(calicost_ok, "ok", "missing_or_empty"),
    stringsAsFactors = FALSE
  )
  ari_rows[[length(ari_rows) + 1]] <- compute_ari_row(
    rep_name = rep_name, tool = "CalicoST",
    gt_df = gt_hard_df, pred_df = calicost_df,
    gt_file = ifelse(gt_hard_ok, gt_hard_out, NA),
    pred_file = ifelse(calicost_ok, calicost_out, NA),
    min_overlap = CFG$ari_min_overlap
  )
}

pred_summary_df <- rbindlist(pred_summary_rows, fill = TRUE)
ari_df <- rbindlist(ari_rows, fill = TRUE)

results_dir <- file.path(opt$root_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

pred_summary_path <- file.path(results_dir, CFG$pred_summary_file)
ari_summary_path <- file.path(results_dir, CFG$ari_summary_file)

fwrite(pred_summary_df, pred_summary_path)
fwrite(ari_df, ari_summary_path)

message(sprintf("\n[Done] Prediction summary: %s", pred_summary_path))
message(sprintf("[Done] ARI summary: %s", ari_summary_path))
