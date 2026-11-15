#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(readr); library(dplyr); library(optparse) })

opt <- parse_args(OptionParser(option_list=list(
  make_option("--histology_csv", type="character"),
  make_option("--epi_labels",   type="character"),
  make_option("--out_keep",     type="character")
)))

detect_delim <- function(path) {
  first <- tryCatch(readr::read_lines(path, n_max = 1), error=function(e) "")
  first <- sub("^\ufeff", "", first)
  if (grepl("\t", first, fixed=TRUE)) "\t" else if (grepl(";", first, fixed=TRUE)) ";" else ","
}
delim <- detect_delim(opt$histology_csv)
h <- suppressMessages(readr::read_delim(opt$histology_csv, delim=delim, col_names=FALSE, show_col_types=FALSE, trim_ws=TRUE))
stopifnot(ncol(h) >= 2)
names(h)[1:2] <- c("barcode","label")
v1 <- tolower(trimws(sub("^\ufeff","", as.character(h$barcode[1]))))
v2 <- tolower(trimws(sub("^\ufeff","", as.character(h$label[1]))))
if (v1 %in% c("barcode","barcodes","spot","cell","cell_id","cellid","feature") ||
    v2 %in% c("label","histology","class","annotation","type")) h <- h[-1, , drop=FALSE]

epi <- unique(trimws(unlist(strsplit(opt$epi_labels, ","))))
keep <- unique(h$barcode[h$label %in% epi])
readr::write_lines(keep, opt$out_keep)
cat(sprintf("[selector] histo-epi keep n=%d -> %s\n", length(keep), opt$out_keep))