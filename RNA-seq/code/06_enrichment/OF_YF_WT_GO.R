
# 1. Load packages --------------------------------------------------------
rm(list=ls());gc()
setwd("/home/adore_org/all_analysis/RNA/GO/")

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)

# 2. Load data ------------------------------------------------------------

# use color
use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')
# input data
data <- fread('./rawdata/WT/OF_YF_WT.csv', sep = ',',header = T, stringsAsFactors = F)


# 3. Analyze --------------------------------------------------------------
diffData <- data


diffData[, type := "ns"]
diffData[log2FoldChange > 0.585 & padj < 0.05, type := "up"][log2FoldChange < -0.585 & padj < 0.05, type := "down"]
table(diffData$type)
  
geneList <- list(
  up = diffData[type == "up", symbol],
  down = diffData[type == "down", symbol]
)
  
egoList <- map(geneList, ~ {
  enrichGO(
    gene = na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = .x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID),
    OrgDb = "org.Mm.eg.db", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)
})

iwalk(egoList, ~ write.csv(.x@result, str_c("./GO/GO_data/",  "OF_YF_WT_FC_0.585_", .y, ".GO.csv")))



