.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# 先装 CRAN 的 scistreer（轻量）
if (!requireNamespace("scistreer", quietly=TRUE)) {
  install.packages("scistreer", repos="https://cloud.r-project.org")
}

# 安装 GitHub 的 hahmmr（核心依赖，不要解析 Suggests）
if (!requireNamespace("hahmmr", quietly=TRUE)) {
  if (!requireNamespace("remotes", quietly=TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org")
  }
  remotes::install_github("kharchenkolab/hahmmr", dependencies=FALSE, upgrade="never")
}

# 安装 numbat 本体（不要解析 Suggests，避免 ggrastr/Cairo/ragg）
if (!requireNamespace("numbat", quietly=TRUE)) {
  remotes::install_github("kharchenkolab/numbat", dependencies=FALSE, upgrade="never")
}

library(numbat); cat("numbat", as.character(packageVersion("numbat")), "\n")

message("[numbat-setup] OK. Installed")