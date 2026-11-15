#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--decon_csv",         type="character"),
  make_option("--epi_celltypes",     type="character"),
  make_option("--epi_threshold",     type="double",  default=0.5),
  make_option("--stroma_celltypes",  type="character"),
  make_option("--imm_celltypes",     type="character"),
  make_option("--ref_top_frac",      type="double",  default=0.10),
  make_option("--out_keep",          type="character"),
  make_option("--out_ref",           type="character")
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

# barcode 列
cand <- tolower(colnames(df))
bidx <- match(TRUE, cand %in% c("barcode","barcodes","cell","spot","cell_id","cellid"))
if (is.na(bidx)) bidx <- 1L
bcol <- colnames(df)[bidx]

sum_by_keys <- function(keys_str){
  keys <- str_split(keys_str, ",")[[1]] |> trimws() |> tolower()
  keys <- keys[nzchar(keys)]
  cols <- colnames(df)[sapply(colnames(df), function(cn){
    cl <- tolower(cn)
    any(str_detect(cl, str_c(keys, collapse="|")))
  })]
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  num <- suppressWarnings(sapply(df[cols], as.numeric))
  if (is.vector(num)) num <- matrix(num, ncol=1)
  rowSums(num, na.rm=TRUE)
}

epi_sum <- sum_by_keys(opt$epi_celltypes)
strom  <- sum_by_keys(opt$stroma_celltypes)
imm    <- sum_by_keys(opt$imm_celltypes)
si_sum <- ifelse(is.na(strom), 0, strom) + ifelse(is.na(imm), 0, imm)

# keep
keep <- df[[bcol]][epi_sum >= opt$epi_threshold] |> as.character() |> unique()

# ref：top X% by si_sum
n_top <- max(1, ceiling(nrow(df) * max(0, min(1, opt$ref_top_frac))))
ord <- order(si_sum, decreasing=TRUE, na.last=NA)
ref <- df[[bcol]][ord][seq_len(min(n_top, length(ord)))] |> as.character() |> unique()

dir.create(dirname(opt$out_keep), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$out_ref),  recursive=TRUE, showWarnings=FALSE)
readr::write_lines(keep, opt$out_keep)
readr::write_lines(ref,  opt$out_ref)

cat(sprintf("[selector] decon keep (epi_sum >= %.3f) n=%d -> %s\n", opt$epi_threshold, length(keep), opt$out_keep))
cat(sprintf("[selector] decon ref  top%.0f%% by stroma+imm n=%d -> %s\n", 100*opt$ref_top_frac, length(ref), opt$out_ref))