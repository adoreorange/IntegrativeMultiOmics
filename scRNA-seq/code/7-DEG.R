# 计算差异表达基因
rm(list=ls());gc()
source('scRNA_scripts/lib.R')
dir.create("6-DEG")
setwd('6-DEG/')
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
sce.all
library(tidyverse)
library(tinyarray)
scRNA = sce.all
head(scRNA@meta.data)
Idents(scRNA) = scRNA$Sample
ct <- c('2MWT', '15MWT', '26MWT')


if (!file.exists('markers_list.Rdata')) {
  markers_list <- list()
  for (i in 1:(length(ct)-1)) {
    for (j in (i+1):length(ct)) {
      markers <- FindMarkers(scRNA, group.by = "Sample",
                             logfc.threshold = 0.1,
                             ident.1 = ct[i],
                             ident.2 = ct[j])
      markers_list[[paste0(ct[i], "_vs_", ct[j])]] <- markers
    }
  }
  save(markers_list,file = 'markers_list.Rdata')
} else {
  
  load('markers_list.Rdata')
}

length(markers_list)
lapply(markers_list,nrow)

all_markers_sig = lapply(markers_list, function(x){
  markers_sig <- subset(x, p_val_adj < 0.01)
})


marker_stat = as.data.frame(lapply(all_markers_sig,function(x){
  # x=all_markers_sig[[1]]
  Up = sum(x$avg_log2FC>0.5)
  Down = sum(x$avg_log2FC< -0.5)
  Total = Up+Down
  return(c(Up, Down, Total))
}))
rownames(marker_stat) = c("Up","Down","Total")
marker_stat
library(gridExtra)

X2MWT_vs_15MWT = all_markers_sig[[1]]
X2MWT_vs_26MWT = all_markers_sig[[2]]
X15MWT_vs_26MWT = all_markers_sig[[3]]
##每一组加上group上下调
deg <- function(dat){
  dat$group = NA
  dat$group[(dat$avg_log2FC < -0.5) & (dat$p_val_adj < 0.01)] <- "down"
  dat$group[(dat$avg_log2FC > 0.5) & (dat$p_val_adj < 0.01)] <- "up"
  table(dat$group )
  dat = na.omit(dat)
  dat
}
P001 = deg(X2MWT_vs_15MWT)
table(P001$group)
P002 = deg(X2MWT_vs_26MWT)
P003 = deg(X15MWT_vs_26MWT)

####veen图####
###上调
X2MWT_vs_15MWT_up = rownames(P001[P001$group == 'up',]) 
X2MWT_vs_26MWT_up = rownames(P002[P002$group == 'up',]) 
X15MWT_vs_26MWT_up = rownames(P003[P003$group == 'up',]) 

inter_upgene <- intersect(intersect(X2MWT_vs_15MWT_up, X2MWT_vs_26MWT_up), X15MWT_vs_26MWT_up)


#三元#
#BiocManager::install("VennDetail")
library(VennDetail)
library(VennDiagram) 
#X2MWT_vs_15MWT,X2MWT_vs_26MWT,X15MWT_vs_26MWT UP
venn <- venndetail(list(X2MWT_vs_15MWT_up = X2MWT_vs_15MWT_up, X2MWT_vs_26MWT_up = X2MWT_vs_26MWT_up, X15MWT_vs_26MWT_up= X15MWT_vs_26MWT_up))
detail(venn) 

# 韦恩图
venn.diagram(x=list(X2MWT_vs_15MWT_up,X2MWT_vs_26MWT_up,X15MWT_vs_26MWT_up),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             category.names = c("IAC_vs_MIA_up", "IAC_vs_AIS_up","MIA_vs_AIS_up") , #标签名
             cat.dist = 0.07, # 标签距离圆圈的远近
             cat.pos = c(-30, -330, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='veenup.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw",# 压缩算法
             height = 2100 ,   # 高度
             width = 2300
             
)
# venn
library(ggvenn)
xx <- list(A = X2MWT_vs_15MWT_up,B = X2MWT_vs_26MWT_up,C=X15MWT_vs_26MWT_up)
p1 <- ggvenn(xx,show_percentage = T,show_elements = F,label_sep = ",",
             digits = 1,stroke_color = "white",
             text_color = "black",text_size = 6,set_name_size = 8,
             fill_color = c("#1E90FF", "#FF8C00","#4DAF4A"),
             set_name_color = c("#1E90FF","#FF8C00","#4DAF4A"))
p1
library(ggpubr)
# 使用Reduce函数来计算列表中所有向量的交集
Intersect <- function(x) {Reduce(intersect, x)}
texttable <- data.frame(`hub-gene` = Intersect(xx))
p2 <- ggtexttable(texttable,rows = NULL)
library(cowplot)
p1 %>% ggdraw() + draw_plot(p2,scale=0.008,x=0.62,y=0.27,width=0.5,height=0.1)
p2

# 韦恩饼图
plot(venn, type = "vennpie")

##多个数据集的韦恩饼图
# vennpie(venn, 
#         min = 4 # 显示集合至少包含来自四个数据集的元素
#         # any = 1, revcolor = "lightgrey" # 突出显示唯一或共享子集
# )

# 韦恩条形图
dplot(venn, order = TRUE, textsize = 4)

# upset图
plot(venn, type = "upset")


# 韦恩饼图
plot(venn, type = "vennpie")

# 韦恩条形图
dplot(venn, order = TRUE, textsize = 4)

# downset图
plot(venn, type = "upset")

##3D PCA
set.seed(98765)
##PCA showing the DEGs among the three groups. The distance between dots represents the difference between groups
table(scRNA$Sample)
meta = scRNA@meta.data
meta$sample_celltype = paste(meta$Sample,meta$celltype,sep = '_')
table(meta$sample_celltype)
#Idents(scRNA)
sce = scRNA
colnames(sce)
rownames(sce)
sce@meta.data = meta
avg <- AverageExpression(object = sce, group.by = "sample_celltype")
a = avg$RNA
b = as.data.frame(a)
###提取出差异基因
DEGgene = c(rownames(X2MWT_vs_15MWT),rownames(X2MWT_vs_26MWT),rownames(X15MWT_vs_26MWT))
##有一些重复的差异基因你
DEG = b[rownames(b) %in% DEGgene,]
str(DEG)
dat = as.data.frame(t(DEG))
pca.res <- prcomp(dat, scale. = T, center = T)
pca.res
tmp <- as.data.frame(pca.res$x)
head(tmp)

library(scatterplot3d)
rownames(tmp)
col = ifelse(str_detect(rownames(tmp),"2MWT"),'purple',
             ifelse(str_detect(rownames(tmp),"15MWT"),'orange','green'))
scatterplot3d(tmp[,1:3],color=col,
              pch = 16,angle=30,
              box=T,type="p",
              lty.hide=2,lty.grid = 2)
legend("topright",c('2MWT','15MWT','26MWT'),
       fill=c('purple','orange','green'),box.col=NA,cex=0.7)




library(ggpubr)
ggscatter(df, x = "MIAAIS", y = "IACMIA",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)

library(ggplot2)
library(ggrepel)
library(ggthemes)
A=df
ggplot(A, aes(x=MIAAIS, y=IACMIA,color = celltype)) +
  geom_hline(yintercept= c(0, 0), color = "black",  size=1) +#添加横线
  geom_vline(xintercept=c(0, 0), color = "black", size=1)+
  geom_point(size = 3,shape=21)

colnames(A)
ggplot(A, aes(x=MIAAIS, y=IACMIA,color = celltype)) +
  geom_hline(yintercept= c(0, 0), color = "black",  size=0.5) +#添加横线
  geom_vline(xintercept=c(0, 0), color = "black", size=0.5)+
  geom_point(size = 3)+
  xlim(-2,3)+
  ylim(-3, 2)+
  labs(x = "Log2FC MIA enriched(MIA vs AIS))",
       y = "Log2FC IAC enriched(IAC vs MIA)", title = "") + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        axis.text = element_text(size = 14, color = 'black'), 
        axis.title = element_text(size = 14, color = 'black'))+
  theme(legend.position = "right")+
  theme_few()+
  geom_text_repel(data=A, aes(label=gene), color="black", size=2, fontface="italic", 
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)


 
###三组差异分析火山图####
library(tidyverse)

all_markers_sig2 = lapply(markers_list, function(x){
  markers_sig <- subset(x,abs(avg_log2FC)>0.5)
})

X2MWT_vs_15MWT = all_markers_sig2[[1]]
X2MWT_vs_26MWT = all_markers_sig2[[2]]
X15MWT_vs_26MWT = all_markers_sig2[[3]]
X2MWT_vs_15MWT$group = 'X2MWT_vs_15MWT'
X2MWT_vs_15MWT$gene = rownames(X2MWT_vs_15MWT)
X2MWT_vs_26MWT$group = 'X2MWT_vs_26MWT'
X2MWT_vs_26MWT$gene = rownames(X2MWT_vs_26MWT)
X15MWT_vs_26MWT$group = 'X15MWT_vs_26MWT'
X15MWT_vs_26MWT$gene = rownames(X15MWT_vs_26MWT)
dat = rbind(X2MWT_vs_15MWT,X2MWT_vs_26MWT,X15MWT_vs_26MWT)

#添加显著性标签：
colnames(dat)
dat$label <- ifelse(dat$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
head(dat)
table(dat$label)
#依次获取最显著的基因
top_15_rows <- dat %>%
  group_by(group) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 15)
table(top_15_rows$group)

#新增一列，将Top差异基因标记为2，其他的标记为1
dat$size <- case_when(!(dat$gene %in% top_15_rows$gene)~ 1,
                      dat$gene %in% top_15_rows$gene ~ 2)
table(dat$size)
#提取非Top10的基因表格；
dt <- filter(dat,size==1)
#绘制散点火山图
dt <- filter(dat,size==1)
head(dt)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = group, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p


#叠加每个Cluster Top10基因散点
str(dt)
dt = as.data.frame(dt)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = group, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top_15_rows,
              aes(x = group, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

#根据图p中log2FC区间确定背景柱长度
dfbar<-data.frame(x=c(1,2,3),
                  y=c(10,10,11))
dfbar1<-data.frame(x=c(1,2,3),
                   y=c(-7.5,-7,-5))
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#把散点火山图叠加到背景柱上
p2 <- ggplot() +
  geom_col(data = dfbar, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_jitter(data = dt, aes(x = group, y = avg_log2FC, color = label), size = 0.85, width = 0.4) +
  geom_jitter(data = top_15_rows, aes(x = group, y = avg_log2FC, color = label), size = 1, width = 0.4) +
  scale_x_discrete() 
p2

#添加X轴的stage色块标签
dfcol<-data.frame(x=c(1:3),
                  y=0,
                  label=c('X15MWT_vs_26MWT','X2MWT_vs_15MWT','X2MWT_vs_26MWT'))
mycol <- c("#00A0877F","#3C54887F","#F39B7F7F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 1.7,
                     show.legend = F)
p3


#给每个stage差异表达前Top15基因加上标签
p4 <- p3+
  geom_text_repel(
    data=top_15_rows,
    aes(x = group, y = avg_log2FC,label = gene),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

#散点颜色调整
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5


#修改X/Y轴标题和添加cluster数字：
p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6

#自定义主题美化：
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7


setwd('../')
