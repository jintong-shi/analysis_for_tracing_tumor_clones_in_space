#!/usr/bin/env Rscript
# Thin R wrapper around extract_groundtruth_subclone.py: resolves a python
# interpreter and runs the carrier-based ground-truth subclone extraction.

suppressPackageStartupMessages({
  library(optparse)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) stop("Cannot resolve script path.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))

option_list <- list(
  make_option("--root_dir", type = "character",
              help = "Root benchmark folder containing Rep_*"),
  make_option("--tool_out_name", type = "character", default = "Tool_Outputs",
              help = "Unused compatibility option [default %default]"),
  make_option("--std_in_name", type = "character", default = "Standardized_Inputs",
              help = "Unused compatibility option [default %default]"),
  make_option("--python_bin", type = "character", default = Sys.getenv("PYTHON_BIN", unset = ""),
              help = "Python executable for the extractor [default: PYTHON_BIN or python3]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$root_dir)) {
  stop("Usage: Rscript extract_groundtruth_subclone.R --root_dir <PATH>")
}
if (!dir.exists(opt$root_dir)) stop("root_dir does not exist: ", opt$root_dir)

conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
python_candidates <- c(
  opt$python_bin,
  if (nzchar(conda_prefix)) file.path(conda_prefix, "envs", "infercnvpy", "bin", "python") else "",
  if (nzchar(conda_prefix)) file.path(dirname(conda_prefix), "infercnvpy", "bin", "python") else "",
  Sys.which("python3"),
  Sys.which("python")
)
python_candidates <- python_candidates[nzchar(python_candidates)]
python_bin <- ""
for (candidate in python_candidates) {
  if (file.exists(candidate) || nzchar(Sys.which(candidate))) {
    python_bin <- candidate
    break
  }
}
if (!nzchar(python_bin)) {
  stop("Cannot find python. Pass --python_bin or set PYTHON_BIN.")
}

extractor <- file.path(script_dir, "extract_groundtruth_subclone.py")
if (!file.exists(extractor)) stop("Missing extractor: ", extractor)

cmd_args <- c(
  extractor,
  "--root-dir", normalizePath(opt$root_dir),
  "--major-frac", "0.05",
  "--carrier-gene-frac", "0.80",
  "--max-cells", "3000"
)

status <- system2(python_bin, cmd_args)
if (!identical(status, 0L)) {
  stop("Ground-truth extractor failed with status: ", status)
}
