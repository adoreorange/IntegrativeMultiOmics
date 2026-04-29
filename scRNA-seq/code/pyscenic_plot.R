## 可视化代码修正版

rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
setwd('/home/hyf/pyscenic/')
source('/home/hyf/pyscenic/code/mycolors.R')

## 提取 out_SCENIC.loom 信息
inputDir='./out_data/'
scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
library(SCENIC)
loom <- open_loom('./out_data/out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

library(SeuratData) # 加载seurat数据集  
sce <- readRDS('/home/hyf/analysis/5-B1/B1_stage.rds')   
table(sce$seurat_clusters)
table(Idents(sce))
Idents(sce) <- sce$RNA_snn_res.0.7
sce$celltype = sce$stage
sce$c <- sce$RNA_snn_res.0.7

#TF 可视化
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)
sce 

#确认是否一致
identical(colnames(sub_regulonAUC), colnames(sce))
#[1] TRUE

# 修正点1：创建正确的cellTypes和cellClusters
cellClusters <- data.frame(row.names = colnames(sce), 
                           seurat_clusters = as.character(sce$seurat_clusters),
                           RNA_snn_res.0.7 = as.character(sce$RNA_snn_res.0.7))

cellTypes <- data.frame(row.names = colnames(sce), 
                        celltype = sce$celltype,
                        RNA_snn_res.0.7 = as.character(sce$RNA_snn_res.0.7))

head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

save(sub_regulonAUC,cellTypes,cellClusters,sce,
     file = 'for_rss_and_visual.Rdata')

####
# 尴尬的是TCF4(+)我这个数据里面没有，换了个PAX5(+)和REL(+)
regulonsToPlot = c('Nfkb2(+)','Myb(+)','Elf4(+)','Junb(+)','Stat1(+)',
                   'Irf1(+)','Irf7(+)','Irf3(+)','Irf2(+)')
regulonsToPlot
sce@meta.data = cbind(sce@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))
Idents(sce) <- sce$celltype
table(Idents(sce))

sce$RNA_snn_res.0.7 <- factor(sce$RNA_snn_res.0.7,levels =  c(0,3,8,7,2,9,5,1,4,6))
Idents(sce) <- sce$RNA_snn_res.0.7
DotPlot(sce, features = unique(regulonsToPlot)) + RotatedAxis()

pdf('AP-1_Ridgeplot_clusters.pdf',height = 8,width = 8)
RidgePlot(sce, features = 'Stat1(+)' , ncol = 1,group.by = 'RNA_snn_res.0.7')
dev.off()

pdf('dotplot_clusters.pdf',height = 8,width = 8)
VlnPlot(sce, features = regulonsToPlot,pt.size = 0,group.by = 'RNA_snn_res.0.7' ) 
dev.off()

VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ) 
FeaturePlot(sce, features = regulonsToPlot )

#### 1. TF活性均值
# Split the cells by cluster:
selectedResolution <- "RNA_snn_res.0.7" # select resolution

# 修正点2：检查列是否存在
if (!selectedResolution %in% colnames(cellClusters)) {
  stop(paste("错误: 列", selectedResolution, "不存在于cellClusters中。可用列:", 
             paste(colnames(cellClusters), collapse = ", ")))
}

# 使用cellClusters进行分组
cellsPerGroup <- split(rownames(cellClusters), 
                       cellClusters[, selectedResolution])

# 检查分组结果
cat("成功创建cellsPerGroup，包含", length(cellsPerGroup), "个组\n")

# 去除extended regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# 检查计算结果
cat("regulonActivity_byGroup维度:", dim(regulonActivity_byGroup), "\n")

# Scale expression. 
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# 2. 热图查看TF分布
# 定义列顺序
col_order <- c('0','3','8','7','2','9','5','1','4','6')

# 确保列顺序正确
if(all(col_order %in% colnames(regulonActivity_byGroup_Scaled))) {
  regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled[, col_order]
} else {
  warning("部分列不存在，使用原始顺序")
}

# 绘制特定IRF家族的热图
pdf('./plot/cluster_order_IRF_Heatmap.pdf', height = 7, width = 9)
pheatmap(regulonActivity_byGroup_Scaled[c('Stat1(+)','Irf1(+)','Irf7(+)','Irf3(+)',"Irf2(+)","Nfkb1(+)"), ], 
         cluster_cols = FALSE,
         fontsize_col = 18, fontsize_row = 16, cellwidth = 50, cellheight = 50)
dev.off()

# 绘制所有regulon的热图
pdf('./plot/cluster_IRF_Heatmap.pdf', height = 25, width = 9)
pheatmap(regulonActivity_byGroup_Scaled, cluster_cols = FALSE,
         fontsize_col = 18, fontsize_row = 16, cellwidth = 10, cellheight = 10)
dev.off()

## 3. rss查看特异TF
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

#### top20
# 计算每行均值
row_means <- rowMeans(regulonActivity_byGroup_Scaled)

# 获取均值最高的20个
top20_idx <- order(row_means, decreasing = TRUE)[1:30]

# 提取并绘图
top20_matrix <- regulonActivity_byGroup_Scaled[top20_idx, ]

pdf('./plot/cluster_order_IRF_top20_mean.pdf', height = 8, width = 8)
pheatmap(top20_matrix, 
         cluster_cols = FALSE,
         fontsize_col = 18, 
         fontsize_row = 10,
         cellwidth = 30, 
         cellheight = 15)
dev.off()

# 4. 其他查看TF方式
library(dplyr) 
rss_matrix <- regulonActivity_byGroup_Scaled
head(rss_matrix)

df <- do.call(rbind,
              lapply(1:ncol(rss_matrix), function(i){
                data.frame(
                  path  = rownames(rss_matrix),
                  cluster = colnames(rss_matrix)[i],
                  sd.1 = rss_matrix[,i],
                  sd.2 = apply(rss_matrix[,-i], 1, median)  
                )
              }))

df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss_matrix[top5$path,] 

annotation_colors <- list(Group = c("2" = "#E5D2DD", "4" = "#53A85F",'7'='#F1BB72','8'='#F3B1A0',
                                    '1'='#D6E7A3','6'='#57C3F3','0'='#476D87','5'='#E95C59','3'='#E59CC4'))

pdf('TF_top5_Heatmap.pdf', height = 12, width = 7)
pheatmap(n,
         annotation_row = rowcn, cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize_col = 14, fontsize_row = 14,
         cellwidth = 15, cellheight = 15, annotation_colors = annotation_colors,
         show_rownames = TRUE)
dev.off()

cat("所有分析完成！\n")