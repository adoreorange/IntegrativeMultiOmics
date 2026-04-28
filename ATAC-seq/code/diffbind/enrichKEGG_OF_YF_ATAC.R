peakAnno_up_df2<-peakAnno_up_df[grepl(x=peakAnno_up_df$annotation,pattern = 'Promoter'),] 
up_gene <- data.frame('gene'=unique(peakAnno_up_df2 $GeneName))
write.csv(up_gene,'up_gene.csv')
ego_KEGG <- enrichKEGG(gene     = peakAnno_up_df2$geneId,
                   #universe     = row.names(dge.celltype),
                   organism      = 'mmu',
                   keyType       = 'kegg',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_kegg <- data.frame(ego_KEGG)
write.csv(ego_kegg,'./ego_kegg_all_OF_YF_WT_UP.csv')


peakAnno_down_df2<-peakAnno_down_df[grepl(x=peakAnno_down_df$annotation,pattern = 'Promoter'),] 
ego_down_KEGG <- enrichKEGG(gene     = peakAnno_down_df2$geneId,
                       #universe     = row.names(dge.celltype),
                       organism      = 'mmu',
                       keyType       = 'kegg',
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
ego_down_kegg <- data.frame(ego_down_KEGG)
write.csv(ego_down,'ego_down.csv')
write.csv(ego_down,'./enrichGO_all_OM_YM_WT_UP.csv')

engo_up_tar <- ego_kegg[which(ego_kegg$ID %in% c('mmu04666','mmu04150','mmu04014',
                                                 'mmu05235','mmu04010','mmu04115',
                                                 'mmu04658','mmu04659','mmu04630',
                                                 'mmu04210','mmu04662','mmu04668')),]

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

