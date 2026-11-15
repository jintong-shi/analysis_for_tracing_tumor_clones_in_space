#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--label_csv", type="character"),
  make_option("--out_ref",   type="character"),
  make_option("--positive",  type="character", default="diploid,normal",
              help="哪些 label 视为正常作为ref，逗号分隔，默认: diploid,normal")
)))

stopifnot(file.exists(opt$label_csv))
pos <- str_split(opt$positive, ",")[[1]] %>% tolower() %>% trimws()

df <- suppressMessages(readr::read_csv(opt$label_csv, show_col_types = FALSE, progress = FALSE))
# 兼容第一列就是barcode的情况
if (!("barcode" %in% names(df))) {
  names(df)[1] <- "barcode"
}
stopifnot("label" %in% names(df))
df <- df %>% mutate(barcode = as.character(barcode),
                    label = tolower(trimws(as.character(label)))) %>%
             filter(label %in% pos) %>%
             distinct(barcode)

dir.create(dirname(opt$out_ref), recursive = TRUE, showWarnings = FALSE)
readr::write_lines(df$barcode, opt$out_ref)
cat(sprintf("[label_to_ref] %s -> %s (n=%d)\n", opt$label_csv, opt$out_ref, nrow(df)))