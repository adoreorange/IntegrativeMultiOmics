# 绘制B-1a细胞样本BCR多样性分析并保存数据
rm(list=ls());gc()
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr) 
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/mycolors.R')
source('./scRNA_scripts/lib.R')
scRNA <- readRDS('./B1/out_data/scRNA_BCR.rds')

combined_sample <- expression2List(sc = scRNA,split.by = 'Sample')

length(combined_sample) #now listed by cluster

p<-quantContig(combined_sample, cloneCall="gene+nt", scale = TRUE,exportTable = T)
write.csv(p,'./sample_plot/sample_quantContig_diversity.csv')
ggsave('./sample_plot/BCR_diversity.pdf',height = 6,width = 7)

compareClonotypes(combined_sample, numbers = 10, 
                  cloneCall="aa", graph = "alluvial")

clonalDiversity(combined_sample, cloneCall = "nt")

p<- clonalHomeostasis(combined_sample, cloneCall = "nt")
write.csv(p,'./sample_plot/sample_clonalHomeostasis.csv')
ggsave('./sample_plot/BCR_clonalHomeostasis.pdf',height = 6,width = 7)


p <-clonalProportion(combined_sample, cloneCall = "nt",exportTable = T)
write.csv(p,'./sample_plot/sample_clonalProportion.csv')
ggsave('./sample_plot/BCR_clonalProportion.pdf',height = 6,width = 7)

quantContig(combined_sample, cloneCall="gene+nt", scale = TRUE) 
ggsave(filename="B1_clonadiveristy_sample.pdf",height = 6,width = 7)

vizGenes(combined_sample,plot ='heatmap', gene = 'J',)
lengthContig(combined_sample, cloneCall="aa", chain = '') 

compareClonotypes(combined_sample, numbers = 10, 
                  cloneCall="aa", graph = "alluvial")



