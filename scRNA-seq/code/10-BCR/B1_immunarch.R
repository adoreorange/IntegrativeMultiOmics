# B-1a细胞BCR分析
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

scRNA <- readRDS('./B1/out_data/scRNA_BCR.rds')
combined_sample <- expression2List(sc = scRNA,split.by = 'Sample')


M2 <- read.csv('./rawdata/WT/2MWT_filtered_contig_annotations.csv',header = T)
M2_sub <- M2[match(combined_sample[['2MWT']]$Barcode,M2$barcode),]
write.csv(M2_sub,'./rawdata/B1_data/2MWT.csv',row.names = F,quote = F)

M15 <- read.csv('./rawdata/WT/15MWT_filtered_contig_annotations.csv',header = T)
M15_sub <- M15[match(combined_sample[['15MWT']]$Barcode,M15$barcode),]
write.csv(M15_sub,'./rawdata/B1_data/15MWT.csv',row.names = F,quote = F)


M26 <- read.csv('./rawdata/WT/26MWT_filtered_contig_annotations.csv',header = T)
M26_sub <- M26[match(combined_sample[['26MWT']]$Barcode,M26$barcode),]
write.csv(M26_sub,'./rawdata/B1_data/26MWT.csv',row.names = F,quote = F)


immdata <- repLoad("./rawdata/B1_data/") 

dir.create('./B1_immu2')
setwd('./B1_immu2/')
immdata$meta$Sample <- factor(immdata$meta$Sample, c('2MWT','15MWT','26MWT'))

exp_vol <- repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c('2MWT','15MWT','26MWT'))

repDiversity(.data = immdata$data, .method = "div", .q = 5, .do.norm = NA, .laplace = 0) %>%
  vis()

imm_pr <- repClonality(immdata$data, .method = "clonal.prop",.perc = 51)
pdf('./B1_clone_prop.pdf',height = 7,width = 8)
vis(imm_pr)
dev.off()

imm_d50 <- repDiversity(immdata$data, .method = 'd50')
pdf('./B1_clone_D50.pdf',height = 7,width = 8)
vis(imm_d50)
dev.off()

imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000))
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000,30000))
#imm_top <- imm_top[match(c('2MWT','15MWT','26MWT'), rownames(imm_top)),]
pdf('./B1_top_clone.pdf',height = 7,width = 8)
vis(imm_top)
dev.off()

imm_rare <- repClonality(immdata$data, .method = "rare")
pdf('./B1_Rare_clone.pdf',height = 7,width = 8)
vis(imm_rare)
dev.off()

imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
pdf('./B1_homeo_clone.pdf',height = 7,width = 8)
vis(imm_hom)
dev.off()

# 各样本的top clone分析
#imm_top <- imm_top[match(c('2MWT','15MWT','26MWT'), rownames(imm_top)),]
vis(imm_top) + vis(imm_top, .by = "Sample", .meta = immdata$meta)
pdf('./B1_top_clone_split.pdf',height = 6,width = 7)
vis(imm_top, .by = "Sample", .meta = immdata$meta)
dev.off()

vis(imm_rare) + vis(imm_rare, .by = "Sample", .meta = immdata$meta)
pdf('./B1_rare_clone_split.pdf',height = 6,width = 7)
vis(imm_rare, .by = "Sample", .meta = immdata$meta)
dev.off()

vis(imm_hom) + vis(imm_hom, .by = 'Sample', .meta = immdata$meta)
pdf('./B1_homeo_clone_split.pdf',height = 6,width = 7)
vis(imm_hom, .by = 'Sample', .meta = immdata$meta)
dev.off()

# 追踪前5个高丰度克隆型
tc <- trackClonotypes(immdata$data, list(1, 10), .col = 'aa+v')
tc$name <- paste0(tc$CDR3.aa, ';', tc$V.name)
tc <- tc [ , c ("name","2MWT" , "15MWT" , "36MWT" ) ]


vis(tc) + 
  scale_fill_manual(values = mycolors) +
  theme_minimal()  # 可选：调整主题

library(immunarch)
library(ggplot2)


# 绘图与美化
set.seed(1)
p <- vis(tc) + 
  scale_fill_manual(
    values = mycolors,
    name = "Clonotypes",          # 修改图例标题
    labels = function(x) gsub("_", " ", x)  # 格式化标签（如去除下划线）
  ) +
  theme_minimal(base_size = 12) +  # 基础字号
  labs(
    title = "Clonotype Frequency Across Samples",
    subtitle = "Target: IGHV11/IGHV12 Clonotypes",
    y = "Frequency (%)"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,          # 标题居中
      margin = margin(b = 10)  # 底部边距
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      color = "grey40"
    ),
    axis.text.x = element_text(
      angle = 45,           # X轴标签倾斜45度
      hjust = 1,            # 水平对齐方式（1=右对齐）
      vjust = 1,            # 垂直对齐方式
      colour = "grey20"
    ),
    axis.text.y = element_text(
      colour = "grey20"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "right",       # 图例置于底部
    legend.title.align = 0.5,          # 图例标题居中
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),  # 去除X轴主网格线
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(          # 控制图例显示方式
                        # 图例单行排列
    title.position = "top",
    keywidth = 1.5,                   # 图例色块宽度
    keyheight = 1.5                   # 图例色块高度
  ));p

# 输出图形
print(p)
ggsave("./B1_immu/clonotype_tracking2.pdf", p, 
       width = 13, height = 6, dpi = 300, bg = "white")
