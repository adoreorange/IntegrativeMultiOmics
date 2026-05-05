# 可视化YM-WT和OM-WT的RRBS数据与RNA-seq数据的四象限DEGs
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
rm(list=ls());gc()
# 1. Read the data
setwd('/home/hyf/rrbs')

data1 <- read.csv("./data/YM-WT_vs_OM-WT-RRBS_bulk-quadrants.csv", row.names = 1)

data <- data1


# 3. set quadrant DEGs

data_significant <- data %>%
  filter(FDR_RRBS < 0.01 & FDR_RNA < 0.01)

quadrant1_genes <- data_significant %>%
  filter(log2FC_RRBS > 10 & log2FC_RNA > 1) %>%
  pull(GeneSymbol)

quadrant2_genes <- data_significant %>%
  filter(log2FC_RRBS < -10 & log2FC_RNA > 1) %>%
  pull(GeneSymbol)

quadrant3_genes <- data_significant %>%
  filter(log2FC_RRBS < -10 & log2FC_RNA < -1) %>%
  pull(GeneSymbol)

quadrant4_genes <- data_significant %>%
  filter(log2FC_RRBS > 10 & log2FC_RNA < -1) %>%
  pull(GeneSymbol)

# 在数据框df中新增quadrant列（需放在绘图代码前执行）
data <- data %>%
  mutate(
    quadrant = case_when(
      GeneSymbol %in% quadrant1_genes ~ "Q1_UpUp",
      GeneSymbol %in% quadrant2_genes  ~ "Q2_DownUp",
      GeneSymbol %in% quadrant3_genes ~ "Q3_DownDown", 
      GeneSymbol %in% quadrant4_genes ~ "Q4_UpDown",
      TRUE ~ "ns"  # 不满足阈值的基因
    )
  )
write.csv(data[(data$quadrant %in% c("Q1_UpUp","Q3_DownDown")), ],'./anaData/YM-WT_vs_OM-WT-RRBS-RNA-quadrants_37.csv')
# 4. select Labeled genes

use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')
up_label <- c('Tnfaip8','Tnfsf4','Map3k15','Cd200','Mndal','Stat1', 'Irf1','Irf3','Irf7','Irf9','Irf8','Ifit1','Ifit2','Ifit3','Oas1','Oasl2','Mx1','Ifitm1','Isg15','Nfkb1','Tgfb1','Ctla4','Il6')
down_label <- c('Dpp4','Slc7a11','Lin28b','Lpcat2','Bmyc','Myc','Ccnd2','Frrs1','Dnaja1','St3gal4','Oxsm','Atp5f1','Sod2','Il1r1','Ackr3')

mycolors <-c('#E64A35','#4DBBD4' ,'#6BD66B','#3C5588'  ,'#F29F80' ,'#01A187',
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
#labelGene <- c(quadrant1_genes,quadrant2_genes,quadrant3_genes,quadrant4_genes)
labelGene <- c(up_label,down_label)
library(ggplot2)
df <- data

###### 
ggplot(df, aes(x = log2FC_RRBS, y = log2FC_RNA))+
  #geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dotted") +
  #geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dotted") +
  geom_point(data = df[df$quadrant %in% "ns", ], 
             shape = 21, color = "grey40", alpha = 0.5, size = 0.2, stroke = 1)+
  geom_point(data = df[!(df$quadrant %in% "ns"), ], aes(x=log2FC_RRBS, y=log2FC_RNA, 
                                                        fill = quadrant, colour = quadrant, 
                                                        size = -log10(FDR_RRBS)), shape = 21, stroke = 0.7) +
  scale_size_continuous(range = c(1, 3)) +
  scale_alpha_continuous(range = c(0.2, 0.85), limits = c(1, 3)) +
  scale_fill_manual(values = c("#E64A35", "#4DBBD4", "#6BD66B", "#F29F80")) +
  scale_color_manual(values = c("#E64A35", "#4DBBD4", "#6BD66B", "#F29F80")) +
  xlim(c(-8, 8)) + ylim(c(-8, 8)) +
  xlab("RRBS_OM_WT_vs_YM_WT(log2FC)")+
  ylab("RNA_OM_WT_vs_YM_WT(log2FC)") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  geom_vline(xintercept = c(-1, 1), color = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "grey80", linewidth = 0.3, linetype = "dashed")

# plot


######
ggplot(df, aes(x = log2FC_RRBS, y = log2FC_RNA))+
  #geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dotted") +
  #geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dotted") +
  geom_point(data = df[df$quadrant %in% "ns", ], 
             shape = 21, color = "grey40", alpha = 0.5, size = 0.2, stroke = 1)+
  geom_point(data = df[!(df$quadrant %in% "ns"), ], aes(x=log2FC_RRBS, y=log2FC_RNA, 
                                                        fill = quadrant, colour = quadrant, 
                                                        size = -log10(FDR_RRBS)), shape = 21, stroke = 0.7) +
  scale_size_continuous(range = c(1, 3)) +
  scale_alpha_continuous(range = c(0.2, 0.85), limits = c(1, 3)) +
  scale_fill_manual(values = c("#E64A35", "#4DBBD4", "#6BD66B", "#F29F80")) +
  scale_color_manual(values = c("#E64A35", "#4DBBD4", "#6BD66B", "#F29F80")) +
  geom_text_repel(
    data= df[df$GeneSymbol %in% labelGene, ],
    aes(x=log2FC_RRBS, y=log2FC_RNA, label = GeneSymbol),
    color = "black",
    size = 4.5,
    box.padding = unit(0.1, "lines"),
    min.segment.length = unit(0.25, "lines"),  # 强制显示所有连线
    max.overlaps = 50,
    segment.color = "black",
    segment.size = 0.3,       # 统一水平偏移
    force = 1,              # 增强排斥力
    max.iter = 5000,        # 提高迭代次数
    seed = 123              # 固定随机种子
  ) +
  xlim(c(-100, 100)) + ylim(c(-10, 10)) +
  xlab("RRBS_B(DNA methylation, YM_WT-OM_WT)")+
  ylab("RNA_OM_WT_vs_YM_WT(log2FC)") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  geom_vline(xintercept = c(-10, 10), color = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "grey80", linewidth = 0.3, linetype = "dashed")
# plot
ggsave('./plot/YM-WT_vs_OM-WT-RRBS_RNA-quadrants.pdf', height = 8, width = 8)

