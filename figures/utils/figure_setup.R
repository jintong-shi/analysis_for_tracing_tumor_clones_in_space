#!/usr/bin/env Rscript
# =====================================================================
# Standalone shared setup for the figures_final Task-1 figures.
#
# This is the ONLY script the task1 figure scripts are allowed to source.
# It provides:
#   - read_csv()           readr-based reader (quiet)
#   - figure_settings() / setting_vector() / setting_order()   settings API
#   - tool_id/tool_label/tool_color + TOOL_TABLE/TOOL_DISPLAY/TOOL_COLORS/
#     TOOL_ALIAS_TABLE      canonical tool naming & colours
#   - dataset_label() + DATASET_LABELS/DATASET_COLORS          canonical datasets
#   - direction_colors/direction_labels                        gain/loss fills
#   - theme_pub()          shared publication theme
#   - skillings_mack()     omnibus test across tools with incomplete blocks
#                          (reused by Fig 2a and Supp S_1A)
#
# It reads figure_settings.tsv from its OWN directory. Unlike the original
# paper/figures/figure_setup.R this carries no task2/task3 project-path globals
# and creates no directories (no side effects), so it is safe to source from a
# fresh figures_final checkout.
# =====================================================================

suppressPackageStartupMessages({
  for (.p in c("dplyr", "tidyr", "purrr", "ggplot2", "scales", "stringr",
               "readr", "tibble")) {
    if (!requireNamespace(.p, quietly = TRUE))
      stop("Required R package is not installed: ", .p, call. = FALSE)
    library(.p, character.only = TRUE)
  }
})

# ---- locate figure_settings.tsv next to this file (works when source()'d) ----
.figure_setup_dir <- (function() {
  for (fr in rev(sys.frames())) {
    if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile, mustWork = TRUE)))
  }
  d <- Sys.getenv("FIG_UTILS_DIR", "")
  if (nzchar(d)) return(normalizePath(d, mustWork = TRUE))
  stop("figure_setup.R: cannot locate its own directory; set FIG_UTILS_DIR.")
})()
FIGURE_SETTINGS_PATH <- file.path(.figure_setup_dir, "figure_settings.tsv")

# ---- bundled font (Arial-metric Liberation Sans) --------------------
# fonts/ next to this file holds LiberationSans*.ttf (Apache-licensed, metric-
# compatible with Arial). Point fontconfig at it via XDG_DATA_HOME so cairo
# devices embed it in the PDF/PNG (a literal "Arial" also resolves to it, as
# Liberation registers as an Arial substitute). Sandbox-safe: fonts travel with
# the repo and the fontconfig cache goes to a writable temp dir. Only set when
# unset, so a host with real system fonts / Arial is left untouched.
.bundled_fonts_dir <- file.path(.figure_setup_dir, "fonts")
if (dir.exists(.bundled_fonts_dir) &&
    length(list.files(.bundled_fonts_dir, pattern = "\\.ttf$", ignore.case = TRUE))) {
  if (Sys.getenv("XDG_DATA_HOME") == "")
    Sys.setenv(XDG_DATA_HOME = .figure_setup_dir)       # fontconfig scans $XDG_DATA_HOME/fonts
  if (Sys.getenv("XDG_CACHE_HOME") == "")
    Sys.setenv(XDG_CACHE_HOME = file.path(tempdir(), "fontconfig_cache"))
  dir.create(Sys.getenv("XDG_CACHE_HOME"), showWarnings = FALSE, recursive = TRUE)
}

read_csv <- function(path, show_col_types = FALSE, ...) {
  readr::read_csv(path, show_col_types = FALSE, ...)
}

.decode_label <- function(x) ifelse(is.na(x), x, gsub("\\\\n", "\n", x))

FIGURE_SETTINGS <- readr::read_tsv(FIGURE_SETTINGS_PATH, na = c("", "NA"),
                                   show_col_types = FALSE) %>%
  mutate(label = .decode_label(.data$label),
         alias = .decode_label(.data$alias))

figure_settings <- function(setting_type = NULL, scope = NULL) {
  out <- FIGURE_SETTINGS
  if (!is.null(setting_type)) out <- out %>% filter(.data$setting_type %in% .env$setting_type)
  if (!is.null(scope))        out <- out %>% filter(.data$scope %in% .env$scope)
  out %>% arrange(.data$order, .data$id)
}

setting_vector <- function(setting_type, scope = "global", value = "label",
                           names_from = "id") {
  rows <- figure_settings(setting_type, scope)
  v <- rows[[value]]; names(v) <- rows[[names_from]]; v
}

setting_order <- function(scope, setting_type = "tool_order", value = "id") {
  figure_settings(setting_type, scope) %>% arrange(.data$order) %>% pull(all_of(value))
}

# ---- tools ----------------------------------------------------------
TOOL_TABLE   <- figure_settings("tool", "global")
TOOL_DISPLAY <- setting_vector("tool", "global", "label")
TOOL_COLORS  <- setting_vector("tool", "global", "color")

TOOL_ALIAS_TABLE <- TOOL_TABLE %>%
  mutate(alias = paste(.data$id, .data$label, .data$alias, sep = "|")) %>%
  separate_rows(alias, sep = "\\|") %>%
  mutate(alias = str_trim(.data$alias)) %>%
  filter(!is.na(.data$alias), .data$alias != "") %>%
  distinct(alias, .keep_all = TRUE)
TOOL_ALIAS_TO_ID <- setNames(TOOL_ALIAS_TABLE$id, TOOL_ALIAS_TABLE$alias)

tool_id <- function(x) {
  o <- unname(TOOL_ALIAS_TO_ID[as.character(x)])
  ifelse(is.na(o), as.character(x), o)
}
tool_label <- function(x) {
  o <- unname(TOOL_DISPLAY[tool_id(x)])
  ifelse(is.na(o), as.character(x), o)
}
tool_color <- function(x) unname(TOOL_COLORS[tool_id(x)])

# ---- datasets + gain/loss direction ---------------------------------
DATASET_LABELS <- setting_vector("dataset", "global", "label")
DATASET_COLORS <- setting_vector("dataset", "global", "color")
dataset_label <- function(x) {
  o <- unname(DATASET_LABELS[as.character(x)])
  ifelse(is.na(o), as.character(x), o)
}
direction_colors <- setting_vector("direction", "synth_task1", "color")
direction_labels <- setting_vector("direction", "synth_task1", "label")

# ---- font resolver --------------------------------------------------
# Return the first family in `prefer` that is actually registered on this
# machine ("sans" always works = device default). Prefer real "Arial"; fall
# back to "Nimbus Sans L", the URW Helvetica/Arial metric-compatible clone
# installed here (systemfonts would otherwise silently substitute DejaVu Sans
# for a literal "Arial", which is NOT Arial). On a host with true Arial this
# returns "Arial" unchanged.
pick_font <- function(prefer = c("Arial", "Liberation Sans", "sans")) {
  fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
  for (f in prefer) if (identical(f, "sans") || f %in% fams) return(f)
  "sans"
}

# ---- shared publication theme (verbatim from paper/figures/figure_setup.R) ----
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size - 1, color = "grey30"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

# ---- shared stats: Skillings-Mack omnibus across tools --------------
# Blocks (reps/samples) x tools, incomplete blocks allowed. Eigen-based
# Moore-Penrose pseudoinverse of the covariance. Identical to the version used
# by the original fig2_compose.R (panel a) and figS_task1A_syn_overall_compose.R.
skillings_mack <- function(value, tool, block) {
  tools  <- sort(unique(as.character(tool)))
  blocks <- sort(unique(as.character(block)))
  M <- matrix(NA_real_, length(blocks), length(tools), dimnames = list(blocks, tools))
  for (i in seq_along(value)) M[as.character(block[i]), as.character(tool[i])] <- value[i]
  M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]
  k <- ncol(M)
  if (k < 2) return(list(statistic = NA_real_, df = NA_integer_, p = NA_real_, k = k))
  A <- setNames(numeric(k), colnames(M))
  S <- matrix(0, k, k, dimnames = list(colnames(M), colnames(M)))
  for (b in seq_len(nrow(M))) {
    row <- M[b, ]; pres <- which(!is.na(row)); kj <- length(pres)
    if (kj < 2) next
    r <- rank(row[pres])
    A[pres] <- A[pres] + sqrt(12 / (kj + 1)) * (r - (kj + 1) / 2)
    for (x in pres) S[x, x] <- S[x, x] + (kj - 1)
    for (x in pres) for (y in pres) if (x != y) S[x, y] <- S[x, y] - 1
  }
  e <- eigen(S, symmetric = TRUE); tol <- max(abs(e$values)) * 1e-8
  dinv <- ifelse(abs(e$values) > tol, 1 / e$values, 0)
  Spi <- e$vectors %*% (dinv * t(e$vectors)); stat <- as.numeric(t(A) %*% Spi %*% A)
  list(statistic = stat, df = k - 1L,
       p = pchisq(stat, k - 1L, lower.tail = FALSE), k = k)
}
