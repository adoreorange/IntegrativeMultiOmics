setwd('/home/adore_org/diffbind/')
OF_YF_peaks_up <- read.csv('./data_out/OF_YF/OF_YF_sig_up.csv')
OF_YF_peaks_up <- peakAnno_all_df

OF_YF_peaks_up_pro<-OF_YF_peaks_up[grepl(x=OF_YF_peaks_up$annotation,pattern ='Promoter'),] 
up_gene <- data.frame('gene'=unique(OF_YF_peaks_up_pro$GeneName))
gene_last <- bitr(geneID = up_gene$gene, 
                  fromType = "SYMBOL",
                  toType = c("ENTREZID", "GENENAME"),
                  OrgDb = 'org.Mm.eg.db')

write.csv(up_gene,'./data_out/OF_YF/OF_YF_up_gene.csv')
ego_UP <- enrichGO(gene          = gene_last$ENTREZID,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'ENTREZID',
                   ont           = "BP",  #设置为ALL时BP, CC, MF都计算
                   pAdjustMethod = "BH", readable=TRUE,
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_up <- data.frame(ego_UP)
write.csv(ego_up,'./data_out/OF_YF/OF_YF_ego_up.csv')

KEGG <- enrichKEGG(gene = gene_last$ENTREZID, 
                   organism = 'mmu', 
                   keyType = 'kegg', 
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.05, 
                   minGSSize = 10)
df_KEGG <- data.frame(KEGG)
dotplot(KEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

engo_up_tar <- ego_up[which(ego_up$Description %in% c('')),]


peakAnno_down_df2<-peakAnno_down_df[grepl(x=peakAnno_down_df$annotation,pattern = 'Promoter'),] 
down_gene <- data.frame('gene'=unique(peakAnno_down_df2$GeneName))
write.csv(down_gene,'down_gene.csv')
ego_DOWN <- enrichGO(gene         = down_gene$gene,
                     #universe     = row.names(dge.celltype),
                     OrgDb         = 'org.Mm.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
ego_down <- data.frame(ego_DOWN)
write.csv(ego_down,'./data_out/OF_YF/ego_down.csv')


engo_up_tar <- ego_up[1:12,]

engo_down_tar <- ego_down[which(ego_down$Description %in% c('regulation of GTPase activity',
                                                            'immune system development','positive regulation of GTPase activity',
                                                            'toll-like receptor 9 signaling pathway','negative regulation of cell migration',
                                                            'negative regulation of phosphorus metabolic process','cell-substrate adhesion',
                                                            'positive regulation of defense response','response to transforming growth factor beta',
                                                            'innate immune response-activating signaling pathway','mmune response-activating signaling pathway',
                                                            'positive regulation of innate immune response','developmental cell growth')),]

# creat up and down data
up.data <- engo_up_tar %>% as.data.frame() %>% select(c('Description','pvalue'))
down.data <- engo_down_tar %>%  as.data.frame() %>% select(c('Description','pvalue'))

p <- go.kegg_plot(up.data,down.data,cols = c('#FFA500','#FF4500'));p
ggsave('./plot/OF_YF/enrich_GO_OF_YF_WT.pdf', p, width = 12, height = 9, dpi = 300)



# plot function ----
go.kegg_plot <- function(up.data,down.data,cols){
  # up.data <- up.data
  # down.data <- down.data
  up.data$group <- 1
  down.data$group <- -1
  dat=rbind(up.data,down.data)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  gk_plot <- ggplot(dat,aes(reorder(Description, pvalue), y=pvalue)) +
    geom_bar(aes(fill=factor((pvalue>0)+1)),stat="identity", width=0.8, position=position_dodge(0.7)) +
    coord_flip() +
    scale_fill_manual(values=cols, guide=FALSE) +
    labs(x="", y="- Log10 (P value)", title = 'GO Term Enrichment') +
    theme_pander()  + ylim(-10,30) +
    theme(panel.grid.major = element_blank(),    # 移除主要网格线
          panel.grid.minor = element_blank(),    # 移除次要网格线
          #axis.ticks.x = element_blank(),       # 用于移除x轴的刻度线
          axis.line.x = element_line(size = 2, colour = "black"), # 设置x轴的线条宽度为2，颜色为黑色
          axis.ticks.length.x = unit(0.20, "cm"), # 设置x轴刻度线的长度为0.20cm
          axis.text.x = element_text(size = 18, margin = margin(t = 0.3, unit = "cm"), hjust = 0), # 设置x轴刻度标签的字体大小为18，上边距为0.3cm，水平对齐方式为0左，0.5中，1右
          axis.ticks.x = element_line(colour = "black",size = 1) , # 设置x轴刻度线的颜色为黑色，大小为1    
          axis.title.x = element_text(size = 16), # 设置x轴标题的字体大小为16
          
          axis.ticks.y = element_blank(),         # 用于移除y轴的刻度线
          axis.text.y  = element_text(size = 16, hjust=1), # 设置y轴刻度标签的字体大小为16，水平对齐方式为1（右对齐）
          axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小为16
          panel.background = element_rect(fill=NULL, colour = 'white'), # 设置图表面板背景为白色，fill=NULL表示背景不填充颜色
          plot.title = element_text(hjust = 0.5, size=16), # 设置图表主标题的水平对齐方式为居中（0.5），字体大小为16
          text = element_text(family = "Times")   # 设置所有文本字体
    )
}

