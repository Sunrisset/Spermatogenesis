# ============== Install packages. ==============

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DropletUtils")

install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran", force = TRUE)
install.packages("ggsci")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("slingshot")

install.packages('Seurat')

options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
install.packages("devtools")

install_github("shenorrLab/cellAlign")

install.packages("imager")
install.packages("scCustomize")
install.packages("paletteer")
install.packages("Hmisc")

install.packages("BiocManager")
BiocManager::install('pcaMethods')
library(devtools)
install_github("velocyto-team/velocyto.R")



install.packages("BASiCS")

install.packages("EnrichedHeatmap")

install.packages("EnsDb.Mmusculus.v79")