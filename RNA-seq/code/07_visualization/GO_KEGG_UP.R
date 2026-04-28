library(Seurat)
library(patchwork)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
rm(list=ls());gc()
setwd('/home/adore_org/Bulk_data_analysis/')

up_gene <- c('Hoxb7','Dlc1','Sdcbp2','Csrnp3','2610037D02Rik','Map3k15','Dmrt2','Mnda','Cyp11a1',
             'Cd59a','C130026I21Rik','Fer','Creb5','Slamf1','Ptpn14','Apbb1','Hus1b','Mid1','Pbx3',
             'Scn4a','Camkk1','Rgs18','Il2rb','Gm1965','Cd3g','Laptm4b','Syt11','Tnfsf4','D10Wsu102e','AI427809',
             'Rcn3','Scimp','Cd72','Serpinb1a','Rasgrp3','Gm6377','Hddc3','Raph1','Tmprss13','Igf1r','Lacc1',
             'Nt5e','Lag3','AW112010','Rom1','Fcgr4','Mcm6','Cd200','Herc3','Aldoc',
             'Tnfaip8','A530032D15Rik','Irf8','Ssbp3','Id3','F11r','Chst3','Arhgef12','Pgs1','Stk3')
features_data <- read.csv('/home/adore_org/Bulk_data_analysis/data_out/OF_YF_WT.csv',header = T, stringsAsFactors = F)
up_features <-features_data$symbol[features_data$change=='up']
down_features <-features_data$symbol[features_data$change=='down']

up_gene<-up_features
print(length(up_gene))
# use color
use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

# enrichGO
ego_UP <- enrichGO(gene          = up_gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_UP_overlap <- data.frame(ego_UP)
ego_UP_overlap <- ego_UP_overlap %>%  mutate(log10pvalue=-log10(pvalue))
write.csv(ego_UP_overlap, './data_out/ego_UP_overlap.csv')

# KEGG 
genelist <- bitr(up_gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist,
                    organism   = 'mmu',
                
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)


ekegg<- setReadable(ekegg,'org.Mm.eg.db','ENTREZID')

write.csv(ekegg@result, './data_out/kegg_UP_overlap.csv')



# plot function ----

gk_plot <- ggplot(head(ego_UP_overlap,10),aes(reorder(Description,log10pvalue), y=log10pvalue)) +
  geom_bar(stat="identity", width=0.8, fill='#FF4500') +
  coord_flip() +
  labs(x="", y="- Log10 (P value)", title = 'OF_OM_overlap_GO') +
  theme_pander()  + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line.x = element_line(size = 2, colour = "black"),#x轴连线
        axis.ticks.length.x = unit(0.20, "cm"),#修改x轴刻度的高度，负号表示向上
        axis.text.x = element_text(size = 18, margin = margin(t = 0.3, unit = "cm"), hjust = 0),##线与数字不要重叠 hjust 0，0.5，1
        axis.ticks.x = element_line(colour = "black",size = 1) ,#修改x轴刻度的线    
        axis.title.x = element_text(size = 16), # 修改x轴文本字体大小
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = 16, hjust=1),
        axis.title.y = element_text(size = 16),
        panel.background = element_rect(fill=NULL, colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=16),
        text = element_text(family = "Times") # 调整字体
  ) 
gk_plot
ggsave('./plot_out/enrich_OF_OM_overlap_GO.pdf', gk_plot, width = 12, height = 9, dpi = 300)

# plot kegg
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count", 
         x="Richfactor", y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}
erich2plot(ekegg_OF_YF@result)
