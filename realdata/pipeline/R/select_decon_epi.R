#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--decon_csv",       type="character"),
  make_option("--epi_celltypes",   type="character"),   # 逗号分隔关键词
  make_option("--epi_threshold",   type="double", default=0.5),
  make_option("--out_keep",        type="character")
)))

strip_bom <- function(x) sub("^\ufeff", "", x)
detect_delim <- function(path){
  first <- tryCatch(readr::read_lines(path, n_max=1), error=function(e) "")
  first <- strip_bom(first)
  if (grepl("\t", first, fixed=TRUE)) "\t" else if (grepl(";", first, fixed=TRUE)) ";" else ","
}

stopifnot(file.exists(opt$decon_csv))
delim <- detect_delim(opt$decon_csv)
df <- suppressMessages(readr::read_delim(opt$decon_csv, delim=delim, show_col_types=FALSE, trim_ws=TRUE))

# 找 barcode 列（宽松）
cand <- tolower(colnames(df))
bidx <- match(TRUE, cand %in% c("barcode","barcodes","cell","spot","cell_id","cellid"))
if (is.na(bidx)) bidx <- 1L
bcol <- colnames(df)[bidx]

# 目标列：部分匹配
keys <- str_split(opt$epi_celltypes, ",")[[1]] |> trimws() |> tolower()
keys <- keys[nzchar(keys)]
stopifnot(length(keys) > 0)

match_cols <- colnames(df)[sapply(colnames(df), function(cn) {
  cl <- tolower(cn)
  any(str_detect(cl, str_c(keys, collapse="|")))
})]
if (length(match_cols) == 0) stop(sprintf("No EPI columns matched by: %s", paste(keys, collapse=",")))

# 数值化并求和
num <- suppressWarnings(sapply(df[match_cols], as.numeric))
if (is.vector(num)) num <- matrix(num, ncol=1)
epi_sum <- rowSums(num, na.rm=TRUE)

keep <- df[[bcol]][epi_sum >= opt$epi_threshold] |> as.character() |> unique()
dir.create(dirname(opt$out_keep), recursive=TRUE, showWarnings=FALSE)
readr::write_lines(keep, opt$out_keep)
cat(sprintf("[selector] decon-epi keep (sum >= %.3f) n=%d -> %s\n", opt$epi_threshold, length(keep), opt$out_keep))