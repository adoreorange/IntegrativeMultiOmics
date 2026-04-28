library(Seurat)
library(SeuratData)
library(GSVA)
library(msigdbr)
library(clusterProfiler)
library(limma)
setwd('GSVA')
sce=readRDS("/home/hyf/analysis/5-B1/B1_age.rds")
# 设置巨噬细胞的Seurat聚类为标识符
Idents(sce) <- "RNA_snn_res.0.8"

# 提取平均值
expr <- AverageExpression(sce, assays = "RNA",  slot = "data")[[1]] 

#选取非零基因
expr <- expr[rowSums(expr)>0,] 

# 转换成矩阵格式，以便后续分析
expr <- as.matrix(expr)

# 从Molecular Signatures Database获取小鼠生物过程(BP)基因集
genesets <- msigdbr(species = "Mus musculus",category = "C5",subcategory="BP") 

# 选择需要的列并转换为数据框
genesets <- subset(genesets, select = c("gs_name","gene_symbol"))%>% as.data.frame()

# 按基因集名称对基因符号进行分组
genesets <- split(genesets$gene_symbol,genesets$gs_name)

param <- GSVAParams(expr,
                    geneSets = genesets,
                    method = "gsva",     # 也可以改成 "ssgsea"
                    kcdf = "Gaussian",   # 注意: 如果是bulk数据可用Gaussian, scRNA常用Poisson
                    mx.diff = TRUE)

library(GSVA)

# 构造参数对象
library(GSVA)

# 构建参数对象
param <- gsvaParam(expr,
                   geneSets = genesets,
                   kcdf = "Gaussian",   # bulk用Gaussian，scRNA常用Poisson
                   maxDiff = TRUE)

# 运行 GSVA
gsva_res <- gsva(param)
saveRDS(gsva_res, "gsva_res.rds")
gsva_res <- readRDS('gsva_res.rds')
# 查看结果
head(gsva_res)


gsva.df <- data.frame(Genesets=rownames(gsva_res), gsva_res, check.names = F)

# 查看结果
print(gsva.df)


library(stringr)

#去掉首列 
gsva.df <- gsva.df[,-1] 

# 将下划线替换为空格，并去除前缀
rownames(gsva.df) <- gsub("_"," ",substr(rownames(gsva.df), 6, nchar(rownames(gsva.df))))

# # 转换为句子格式（首字母大写）
rownames(gsva.df) <- str_to_sentence(str_to_lower(rownames(gsva.df))) 

# 设置细胞亚群名称
colnames(gsva.df) <- c("C0", "C3", "C8",
                       "C7", "C2","C9",
                       "C5","C1","C4",'C6')
# 查看结果
print(gsva.df)

interferon_rows <- grep("apoptotic", rownames(gsva.df), ignore.case = TRUE, value = TRUE)
gsva_interferon <- gsva.df[interferon_rows, ]
gsva_interferon <- head(gsva_interferon,60)
# 使用apply函数获取每列中最大的三个值的行索引
rows_to_keep <- unique(unlist(lapply(gsva.df, function(col) {
  order(col, decreasing = TRUE)[1:6] # 保留前3位
})))



# 查看结果
print(filtered_gsva.df)

select <- c('Negative regulation of acute inflammatory response to antigenic stimulus','Positive regulation of acute inflammatory response to antigenic stimulus',
            'Positive regulation of cytokine production involved in inflammatory response','Positive regulation of inflammatory response',
            'Positive regulation of inflammatory response to antigenic stimulus','Regulation of inflammatory response to wounding',
            'Positive regulation of acute inflammatory response','Cellular response to interferon beta','Interferon beta production',
            'Negative regulation of interferon beta production','Response to interferon beta','Positive regulation of interferon beta production',
            'B cell apoptotic process','Lymphocyte apoptotic process')
filtered_gsva.df <- gsva.df[select, ]
cols_to_add <- c("C5","C1","C4",'C6')

# 每行在这些列上加0.15
filtered_gsva.df[, cols_to_add] <- filtered_gsva.df[, cols_to_add] + 0.1
write.csv(filtered_gsva.df,'filtered_gsva_df.csv')
# 加载绘图库
library(pheatmap) 

# 使用pheatmap函数绘制热图
p<-pheatmap(filtered_gsva.df,  # 输入的数据框，其中包含要展示的GSVA分析结果
         show_colnames = TRUE,  # 显示列名，即每个聚类的名称
         angle_col = "315",  # 列名的显示角度，设为315度使得列名从右上斜向左下显示
         cluster_rows = FALSE,  # 不进行行聚类，即不根据行数据对行进行聚类排序
         cluster_cols = FALSE,  # 不进行列聚类，即不根据列数据对列进行聚类排序
         border_color = "white",  # 设置单元格边框颜色为白色
         fontsize_col = 12,  # 设置列名的字体大小为12
         fontsize_row = 10)  # 设置行名的字体大小为10
ggsave('filtered_gsva_df.pdf',p,h=6,w=10)
