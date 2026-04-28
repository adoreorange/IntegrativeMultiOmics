library(Seurat)
library(patchwork)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
rm(list=ls());gc()
setwd('/home/adore_org/Bulk_data_analysis/')

down_gene <- c('Adm','Emb','Frrs1','Dnaja1','Dpp4','Hbb-b1','Pde2a','BC035044',
               '1700017B05Rik','Cd209a','Amz1','Eps8','Edaradd','Hsph1',
               'Ehd3','Znrf3','Ctnna1','Alpl','Gpr137b','Tcf7','Samd8','Retnla','Gpr171',
               'Fam78a','St3gal4','Lrrc16a','Serpine2','Scamp1',
               'Cecr2','Otub2','Dag1','Oxsm','Plxnb2','Cd93','Wee1','Rasgef1a',
               'Nrep','2810025M15Rik','Fam134b','Gpr132',
               'Rps4l','Gsto1','Mgst2','Slc16a1','Gprc5b','Il9r')
features_data <- read.csv('/home/adore_org/Bulk_data_analysis/data_out/OM_YM_WT.csv',header = T, stringsAsFactors = F)
up_features <-features_data$symbol[features_data$change=='up']
down_features <-features_data$symbol[features_data$change=='down']

down_gene<-down_features
print(length(down_gene))

# use color
use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

# enrichGO
ego_DOWN <- enrichGO(gene       = down_gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ego_DOWN_overlap <- data.frame(ego_DOWN)
ego_DOWN_overlap <- ego_DOWN_overlap %>%  mutate(log10pvalue=-log10(pvalue))
write.csv(ego_DOWN_overlap, './data_out/ego_DOWN_overlap.csv')

# KEGG 
genelist <- bitr(down_gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')

genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist,
                    organism   = 'mmu',
                    keyType = "kegg",
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
