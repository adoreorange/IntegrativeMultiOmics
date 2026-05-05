# 绘制B-2细胞样本BCR多样性分析并保存数据
rm(list=ls());gc()
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr) 
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/mycolors.R')
source('./scRNA_scripts/lib.R')

setwd("/home/adore_org/B_scRNA-seq/analysis/BCR_data/B2/")

combined <- readRDS('/home/hyf/analysis/BCR_data/B1/out_data/combined-3WT.rds')
sce <-readRDS('/home/adore_org/B_scRNA-seq/analysis/6-B2/2-harmony_B2/sce.all_int.rds')

combined_sample <- expression2List(sc = scRNA,split.by = 'Sample')

length(combined_sample) #now listed by cluster

p<-quantContig(combined_sample, cloneCall="gene+nt", scale = TRUE)
write.csv(p,'./sample_plot/sample_quantContig_diversity.csv')
ggsave('./sample_plot/BCR_diversity.pdf',height = 6,width = 7)



clonalDiversity(combined_sample, cloneCall = "nt")

p<- clonalHomeostasis(combined_sample, cloneCall = "nt")
write.csv(p,'./sample_plot/sample_clonalHomeostasis.csv')
ggsave('./sample_plot/BCR_clonalHomeostasis.pdf',height = 6,width = 7)


p <-clonalProportion(combined_sample, cloneCall = "nt")
write.csv(p,'./sample_plot/sample_clonalProportion.csv')
ggsave('./sample_plot/BCR_clonalProportion.pdf',height = 6,width = 7)

quantContig(combined_sample, cloneCall="gene+nt", scale = TRUE) 
ggsave(filename="B1_clonadiveristy_sample.pdf",height = 6,width = 7)

vizGenes(combined_sample,plot ='heatmap', gene = 'J',)
lengthContig(combined_sample, cloneCall="aa", chain = '') 

compareClonotypes(combined_sample, numbers = 10, 
                  cloneCall="aa", graph = "alluvial")



