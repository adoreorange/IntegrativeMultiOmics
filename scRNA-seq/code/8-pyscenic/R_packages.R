# 检查并安装所需包----
packages <- c("GenomicFeatures","AUCell", "RcisTarget", "GENIE3", "zoo", "mixtools", "rbokeh",
              "DT", "NMF", "pheatmap", "R2HTML", "Rtsne", "doMC", "doRNG")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# 对于从github安装的包，使用devtools ----
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("aertslab/SCopeLoomR")
if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
  devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
}
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

# 安装senic
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC")
packageVersion("SCENIC")
## SCENIC需要一些依赖包，先安装好
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3"))
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))
BiocManager::install(c("doMC", "doRNG"))
BiocManager::install(c("optparse"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(optparse)
