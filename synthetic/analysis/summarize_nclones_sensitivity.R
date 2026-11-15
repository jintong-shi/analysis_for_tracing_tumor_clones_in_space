#!/usr/bin/env Rscript
# Aggregate STARCH + CalicoST subclone ARI across n_clones and compare to n=4.
# Writes subclone_ari_{by_nclones,nclones_summary,delta_vs_n4}.csv to <sweep_dir>/results/.
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sweep_dir",    type = "character"),
  make_option("--baseline_ari", type = "character", default = ""),
  make_option("--baseline_n",   type = "integer",   default = 4L)
)))
stopifnot(!is.null(opt$sweep_dir))
TOOLS <- c("STARCH", "CalicoST")

read_one <- function(path, n) {
  if (!file.exists(path)) { message("skip (missing): ", path); return(NULL) }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$tool %in% TOOLS, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$n_clones_param <- n
  df
}

# per-K sweep outputs (n=4 comes from baseline, not the sweep roots)
parts <- list()
for (K in c(2, 3, 5, 6)) {
  parts[[length(parts) + 1]] <-
    read_one(file.path(opt$sweep_dir, "roots", sprintf("nclones_%d", K),
                       "results", "subclone_ari_summary.csv"), K)
}
if (nzchar(opt$baseline_ari))
  parts[[length(parts) + 1]] <- read_one(opt$baseline_ari, opt$baseline_n)

parts <- Filter(Negate(is.null), parts)
if (length(parts) == 0) stop("No ARI summaries found — nothing to aggregate.")
cols <- Reduce(intersect, lapply(parts, names))
by_n <- do.call(rbind, lapply(parts, function(d) d[, cols, drop = FALSE]))
by_n <- by_n[order(by_n$tool, by_n$n_clones_param, by_n$rep), ]

res_dir <- file.path(opt$sweep_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
write.csv(by_n, file.path(res_dir, "subclone_ari_by_nclones.csv"), row.names = FALSE)

# ---- per tool x n_clones summary: ARI (over scored reps) + monoclone accuracy
grp <- split(by_n, list(by_n$tool, by_n$n_clones_param), drop = TRUE)
summ <- do.call(rbind, lapply(grp, function(g) {
  ari  <- suppressWarnings(as.numeric(g$ari))
  ok   <- !is.na(ari)
  mono <- g$n_gt_subclones == 1                      # ground-truth monoclone reps
  npred <- suppressWarnings(as.numeric(g$n_pred_subclones))
  mono_correct <- mono & !is.na(npred) & npred == 1  # tool also called a single clone
  data.frame(
    tool            = g$tool[1],
    n_clones_param  = g$n_clones_param[1],
    n_reps          = nrow(g),
    n_scored        = sum(ok),
    mean_ari        = if (any(ok)) mean(ari[ok]) else NA_real_,
    median_ari      = if (any(ok)) median(ari[ok]) else NA_real_,
    n_monoclone_gt  = sum(mono),
    monoclone_acc   = if (sum(mono) > 0) sum(mono_correct) / sum(mono) else NA_real_,
    stringsAsFactors = FALSE)
}))
summ <- summ[order(summ$tool, summ$n_clones_param), ]
write.csv(summ, file.path(res_dir, "subclone_ari_nclones_summary.csv"), row.names = FALSE)

# ---- paired delta vs n=4 (per rep x tool) -----------------------------------
n4 <- by_n[by_n$n_clones_param == opt$baseline_n,
           c("rep", "tool", "ari", "n_pred_subclones")]
names(n4) <- c("rep", "tool", "ari_n4", "n_pred_subclones_n4")
delta <- merge(by_n[, c("rep", "tool", "n_clones_param", "ari",
                        "n_pred_subclones", "n_gt_subclones", "status")],
               n4, by = c("rep", "tool"), all.x = TRUE)
delta$ari_delta_vs_n4 <-
  suppressWarnings(as.numeric(delta$ari) - as.numeric(delta$ari_n4))
delta <- delta[order(delta$tool, delta$n_clones_param, delta$rep), ]
write.csv(delta, file.path(res_dir, "subclone_ari_delta_vs_n4.csv"), row.names = FALSE)

cat("\n=== subclone_ari_nclones_summary ===\n")
print(summ, row.names = FALSE)
cat("\nWrote 3 CSVs to:", res_dir, "\n")
