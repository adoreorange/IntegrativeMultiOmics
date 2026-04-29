
library(patchwork) # 拼图
library(ggplot2) # 绘图
library(Seurat) # readRDS, VlnPlot
library(scales) # show_col
library(tidyr)
colors <- c("#0072B2","#009E73","#D55E00","#CC79A7","#F0E442",
            "#56B4E9","#E69F00","#00ADA9","#D0E429","#ED008C","#68217A","#739B57", "#EFE685",
            "#446983", "#BB6239", "#4F4FFF", "#7F2268", "#800202", "#D8D8CD")
my_stacked_violin_plot=function(gene,seurat_object,Clusters='seurat_clusters',
                             width=13,height=10.3,filename="",text.size=10,
                             col=colors){
  sample_sce <- seurat_object[gene,]
  data <- data.frame(t(sample_sce@assays$RNA@data))
  data$Barcode <- rownames(data)
  data$Idents <- sample_sce@meta.data[,Clusters]
  data2 <- gather(data = data,key = 'Gene',value = 'Expr',-c(Barcode,Idents))
  data2$Gene <- factor(data2$Gene, levels = levels(gene))
  p <- ggplot(data2, aes(factor(Idents), Expr, fill = Idents)) + # 基础映射
    # 小提琴图, 每一个参数都有意义
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    # y 轴设置
    scale_y_continuous(expand = c(0, 0),
                       position="right",
                       labels = function(x)
                         c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) + 
    # 分面
    facet_grid(rows = vars((Gene)), scales = "free", switch = "y") +
    # 自定义填充色
    scale_fill_manual(values = colors) +
    # 使用 cowplot 主题
    theme_cowplot(font_size = 12) + 
    # 进一步设置主题，特别是与分面相关的参数要注意哈
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0),
          axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0.5, size = 8)) +
    xlab("Gene Features") + # x 轴表题
    ylab("Expression Level") # y 轴表题
  pdf(paste0(filename,"_vlnplot",".pdf"),width = width,height = height)
  print(p)
  dev.off()
}






