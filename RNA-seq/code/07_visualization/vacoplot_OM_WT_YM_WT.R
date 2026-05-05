library(Seurat)
library(patchwork)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
rm(list=ls());gc()
setwd('/home/adore_org/Bulk_data_analysis/')

# use color
use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')
# input data
data <- fread('/home/adore_org/Bulk_data_analysis/DEGs/OM-WT_YM-WT.csv', sep = ',',header = T, stringsAsFactors = F)

deg_all <- dplyr::select(data,c('V1','log2FoldChange_DESeq2','padj_DESeq2'))
deg_all <- as.data.frame(deg_all)
rownames(deg_all) <- deg_all$V1
deg_all <- deg_all[-1]
# creat dataframe
dif2=data.frame(
  symbol=rownames(deg_all),
  log2FoldChange=deg_all$log2FoldChange_DESeq2,
  padj=deg_all$padj_DESeq2
)

# analysis
logFC = 1
P.Value = 0.05
k1 <- (dif2$padj < P.Value) & (dif2$log2FoldChange < -logFC)
k2 <- (dif2$padj < P.Value) & (dif2$log2FoldChange > logFC)
dif2 <- mutate(dif2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(dif2$change)

# write xlsx
dif2 <- dif2 %>%  mutate(log10padj=-log10(padj))
# write.xlsx(dif2,'./data_out/OF_YF_WT.xlsx')
write.csv(dif2,'./data_out/OM_YM_WT.csv')
dif2 <- dif2[]
up <- length(which(k1=='TRUE'))
down <- length(which(k2=='TRUE'))
title <- paste0('pvalue:', P.Value,';log2FC:',logFC,';Up:',up,';Down:',down,';Total:',length(dif2$symbol))
# select gene
genes <- c()

# plot
p <- VolcanoPlot(dif2, padj=0.05, title=title, label.max = 40, label.symbols = genes,log2FC = 1,cols=c('#FFA500', '#FF4500'))
ggsave('./plot_out/Volcano_OM_YM_WT.pdf', p, width = 10, height = 10, dpi = 600)

# plot function ----
VolcanoPlot=function(dif, log2FC=log2(1.5), padj=0.05, 
                     label.symbols=NULL, label.max=30, 
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="stable";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','stable','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=1, size=3) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color=cols)+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+ ylim(-1, 52) +
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10),
          text = element_text(family = "Times")
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'stable',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "stable" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print((dd_text))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1,
                       segment.color = "black",
                       segment.size = 0.5,      # 设置连线粗细
                       segment.alpha = 0.5,     # 设置连线透明度
                       box.padding = 0.5,       # 设置文本框的内边距
                       force = 1.25,               # 设置排斥力
                       max.overlaps = Inf
  )
}
