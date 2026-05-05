# B-2细胞BCR分析
# 运行immunarch分析
# @param scRNA scRNA对象
# @return immdata对象
setwd('/home/hyf/analysis/BCR_data/')
rm(list=ls());gc()
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr) 
library(immunarch)
library(scRepertoire)

scRNA <- readRDS('./B2/out_data/B2_BCR_1205.rds')
combined_sample <- expression2List(sc = scRNA,split.by = 'Sample')



M2 <- read.csv('./rawdata/WT/2MWT_filtered_contig_annotations.csv',header = T)
M2_sub <- M2[match(combined_sample[['2MWT']]$Barcode,M2$barcode),]
write.csv(M2_sub,'./B2_immu_1025/data/2MWT.csv',row.names = F,quote = F)

M15 <- read.csv('./rawdata/WT/15MWT_filtered_contig_annotations.csv',header = T)
M15_sub <- M15[match(combined_sample[['15MWT']]$Barcode,M15$barcode),]
write.csv(M15_sub,'./B2_immu_1025/data/15MWT.csv',row.names = F,quote = F)


M26 <- read.csv('./rawdata/WT/26MWT_filtered_contig_annotations.csv',header = T)
M26_sub <- M26[match(combined_sample[['26MWT']]$Barcode,M26$barcode),]
write.csv(M26_sub,'./B2_immu_1025/data/26MWT.csv',row.names = F,quote = F)


immdata <- repLoad("./B2_immu2/data/")
immdata <- repLoad("./rawdata/B2_data/") 
dir.create('./B2_immu2')
setwd('./B2_immu2/')

immdata$meta$Sample <- factor(immdata$meta$Sample, c('2MWT','15MWT','26MWT'))

exp_vol <- repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c('2MWT','15MWT','26MWT'))

repDiversity(.data = immdata$data, .method = "div", .q = 5, .do.norm = NA, .laplace = 0) %>%
  vis()

imm_pr <- repClonality(immdata$data, .method = "clonal.prop",.perc = 50)
pdf('./B2_clone_prop.pdf',height = 7,width = 8)
vis(imm_pr,.by = 'Sample',.meta = immdata$meta)
dev.off()

imm_d50 <- repDiversity(immdata$data, .method = 'd50',.perc = 50)
pdf('./B2_clone_D50.pdf',height = 7,width = 8)
vis(imm_d50,.meta = immdata$meta)
dev.off()


imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000),)
#imm_top <- imm_top[match(c('2MWT','15MWT','26MWT'), rownames(imm_top)),]
pdf('./B2_top_clone.pdf',height = 7,width = 8)
vis(imm_top)
dev.off()

imm_rare <- repClonality(immdata$data, .method = "rare")
pdf('./B2_Rare_clone.pdf',height = 7,width = 8)
vis(imm_rare)
dev.off()

imm_hom <- repClonality(immdata$data, 
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
pdf('./B2_homeo_clone.pdf',height = 7,width = 8)
vis(imm_hom)
dev.off()

# 各样本的top clone分析
#imm_top <- imm_top[match(c('2MWT','15MWT','26MWT'), rownames(imm_top)),]
vis(imm_top) + vis(imm_top, .by = "Sample", .meta = immdata$meta)
pdf('./B2_top_clone_split.pdf',height = 6,width = 7)
vis(imm_top, .by = "Sample", .meta = immdata$meta)
dev.off()

vis(imm_rare) + vis(imm_rare, .by = "Sample", .meta = immdata$meta)
pdf('./B2_rare_clone_split.pdf',height = 6,width = 7)
vis(imm_rare, .by = "Sample", .meta = immdata$meta)
dev.off()

vis(imm_hom) + vis(imm_hom, .by = 'Sample', .meta = immdata$meta)
pdf('./B2_homeo_clone_split.pdf',height = 6,width = 7)
vis(imm_hom, .by = 'Sample', .meta = immdata$meta)
dev.off()
