#!/usr/bin/env Rscript
pkgs <- c("remotes","BiocManager","keras","tensorflow")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
msg <- function(...) cat("[r-tuscan setup]", ..., "\n")

# 主包
if (!requireNamespace("TUSCAN", quietly=TRUE)) {
  msg("Installing TUSCAN from GitHub CZang409/TUSCAN ...")
  remotes::install_github("CZang409/TUSCAN", upgrade = "never")
}

msg("Installing Keras (python backend) via keras3 ...")
library(keras)
install_keras()

msg("r-tuscan packages ready.")