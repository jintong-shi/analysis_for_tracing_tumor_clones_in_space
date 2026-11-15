#!/usr/bin/env Rscript
pkgs <- c("remotes","BiocManager")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")

msg <- function(...) cat("[r-shared setup]", ..., "\n")

# infercnv 已由 conda 安装；这里只处理 GitHub 包
need_install <- function(pkg) !requireNamespace(pkg, quietly=TRUE)

if (need_install("yaGST")) {
  msg("Installing yaGST from GitHub miccec/yaGST ...")
  remotes::install_github("miccec/yaGST", upgrade = "never")
}

if (need_install("copykat")) {
  msg("Installing copykat from GitHub navinlabcode/copykat ...")
  remotes::install_github("navinlabcode/copykat", upgrade = "never")
}

if (need_install("SCEVAN")) {
  msg("Installing SCEVAN from GitHub AntonioDeFalco/SCEVAN ...")
  remotes::install_github("AntonioDeFalco/SCEVAN", upgrade = "never")
}

# --- CONICSmat deps & package ---
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  pkgs_cran <- c("beanplot","mixtools","pheatmap","zoo","squash","remotes")
  for (p in pkgs_cran) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt", update = FALSE, ask = FALSE)
  }
  # CONICSmat from GitHub
  if (!requireNamespace("CONICSmat", quietly = TRUE)) {
    remotes::install_github("diazlab/CONICS/CONICSmat", dep = FALSE, upgrade = "never")
  }
})

msg("r-shared packages ready.")