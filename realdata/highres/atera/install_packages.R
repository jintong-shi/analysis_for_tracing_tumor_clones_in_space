dir.create("~/Rlibs/xenium_R442", recursive = TRUE, showWarnings = FALSE)
.libPaths(c("~/Rlibs/xenium_R442", .libPaths()))

install.packages(c("remotes", "BiocManager"))

# 让 Bioconductor 进入和 R 4.4 匹配的 release
BiocManager::install(version = "3.20", ask = FALSE)

# CRAN / archive 指定版本
remotes::install_version("SeuratObject", version = "5.1.0", upgrade = "never")
remotes::install_version("Seurat",       version = "5.3.1", upgrade = "never")
remotes::install_version("ggplot2",      version = "3.5.2", upgrade = "never")
remotes::install_version("future",       version = "1.67.0", upgrade = "never")
remotes::install_version("sp",           version = "2.2-0", upgrade = "never")

# Bioconductor 包指定安装
BiocManager::install(c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "BiocGenerics",
  "GenomeInfoDb",
  "SummarizedExperiment",
  "DelayedArray",
  "Biobase",
  "SparseArray",
  "MatrixGenerics",
  "XVector",
  "zlibbioc"
), ask = FALSE, update = FALSE)