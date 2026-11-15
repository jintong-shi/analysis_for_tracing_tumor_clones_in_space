#!/usr/bin/env Rscript

# 用法:
# Rscript R/select_histology_labels.R \
#   --histology_csv <path> \
#   --labels "Benign,Cancer 3+4" \
#   --out_keep <path/to/keep_barcodes.txt>
#
# 说明:
# - 自动猜测条形码列(barcode列)和标签列(label列)。
# - 条形码列优先匹配: barcode, barcodes, spot, spot_id, barcode_id (忽略大小写)。
# - 标签列会优先选择与 --labels 值命中最多的字符列；
#   若没有任何命中，再按常见列名匹配 label/histology/final_annotations 等。

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0) return(NULL)
  if (i == length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

hist_csv <- get_arg("--histology_csv")
labels_str <- get_arg("--labels")
out_keep <- get_arg("--out_keep")

if (is.null(hist_csv) || is.null(labels_str) || is.null(out_keep)) {
  cat("ERROR: missing required args.\n", file=stderr())
  quit(status = 2)
}

# 读表
df <- tryCatch(
  {
    read.csv(hist_csv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  },
  error = function(e) {
    cat("ERROR: failed to read histology_csv:", conditionMessage(e), "\n", file=stderr())
    quit(status = 2)
  }
)

# 拆分标签
target_labels <- trimws(strsplit(labels_str, ",")[[1]])
target_lc <- tolower(target_labels)

# 猜 barcode 列
guess_barcode_col <- function(nms) {
  prefs <- c("barcode","barcodes","spot","spot_id","barcode_id","cell","cell_id")
  for (p in prefs) {
    hit <- which(tolower(nms) == tolower(p))
    if (length(hit)) return(nms[hit[1]])
  }
  # fallback: 第一个字符型且唯一值数较大的列
  char_cols <- nms[sapply(df, is.character)]
  if (length(char_cols) == 0) return(nms[1])
  uniq_cnt <- sapply(df[char_cols], function(x) length(unique(x)))
  char_cols[which.max(uniq_cnt)]
}

# 猜 label 列
guess_label_col <- function(nms) {
  char_cols <- nms[sapply(df, is.character)]
  if (length(char_cols) == 0) return(nms[1])

  scores <- sapply(char_cols, function(col) {
    sum(tolower(trimws(df[[col]])) %in% target_lc, na.rm = TRUE)
  })
  if (max(scores) > 0) {
    tied <- char_cols[scores == max(scores)]
    prefs <- c("label","labels","final_annotations","final_annotation","final_label",
               "annotation","histology","histo","class","type")
    for (p in prefs) {
      hit <- which(tolower(tied) == tolower(p))
      if (length(hit)) return(tied[hit[1]])
    }
    return(tied[1])
  }

  prefs <- c("label","labels","final_annotations","final_annotation","final_label",
             "annotation","histology","histo","class","type")
  for (p in prefs) {
    hit <- which(tolower(nms) == tolower(p))
    if (length(hit)) return(nms[hit[1]])
  }
  char_cols[1]
}

nms <- names(df)
barcode_col <- guess_barcode_col(nms)
label_col   <- guess_label_col(nms)

if (!(barcode_col %in% nms) || !(label_col %in% nms)) {
  cat("ERROR: failed to guess barcode/label column.\n", file=stderr())
  quit(status = 2)
}

sel <- which(tolower(trimws(df[[label_col]])) %in% target_lc)
barcodes <- unique(trimws(df[[barcode_col]][sel]))
barcodes <- barcodes[nzchar(barcodes)]

# 写出
odir <- dirname(out_keep)
if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
con <- file(out_keep, open = "w", encoding = "UTF-8")
writeLines(barcodes, con = con, sep = "\n", useBytes = TRUE)
close(con)

cat(sprintf("[select_histology_labels] histology_csv=%s | labels=%s | barcode_col=%s | label_col=%s | kept=%d\n",
            hist_csv, labels_str, barcode_col, label_col, length(barcodes)))
