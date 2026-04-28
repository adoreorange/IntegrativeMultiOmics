#!/usr/bin/env Rscript
# DESeq2 差异表达分析
# 用法: Rscript deseq2_analysis.R <count_matrix> <sample_info> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
count_file <- ifelse(length(args) >= 1, args[1], "output/04_quantification/count_matrix.txt")
sample_info_file <- ifelse(length(args) >= 2, args[2], "data/sample_info.txt")
output_dir <- ifelse(length(args) >= 3, args[3], "output/05_de_analysis")

cat("==========================================\n")
cat("DESeq2 差异表达分析\n")
cat("==========================================\n")
cat("计数矩阵:", count_file, "\n")
cat("样本信息:", sample_info_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取数据
cat("读取计数矩阵...\n")
count_matrix <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
count_matrix <- as.matrix(count_matrix)

# 过滤低表达基因
cat("过滤低表达基因...\n")
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, ]
cat("保留基因数:", nrow(count_matrix), "\n")

# 读取样本信息
cat("读取样本信息...\n")
sample_info <- read.table(sample_info_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# 确保样本信息与计数矩阵匹配
common_samples <- intersect(colnames(count_matrix), rownames(sample_info))
count_matrix <- count_matrix[, common_samples]
sample_info <- sample_info[common_samples, , drop = FALSE]

cat("样本数量:", ncol(count_matrix), "\n")
cat("样本分组:", unique(sample_info$condition), "\n")

# 创建 DESeq2 数据集
cat("创建 DESeq2 数据集...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# 运行 DESeq2
cat("运行 DESeq2 分析...\n")
dds <- DESeq(dds)

# 获取结果
cat("提取差异表达结果...\n")
res <- results(dds, alpha = 0.05)
res <- lfcShrink(dds, coef = 2, type = "apeglm")

# 排序并保存结果
res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), file = file.path(output_dir, "DEG_results_all.csv"))

# 筛选显著差异基因
deg <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(deg), file = file.path(output_dir, "DEG_significant.csv"))

cat("总基因数:", nrow(res), "\n")
cat("显著差异基因数:", nrow(deg), "\n")

# ==================== 可视化 ====================

# MA 图
cat("生成 MA 图...\n")
png(file.path(output_dir, "MA_plot.png"), width = 800, height = 600)
plotMA(res, main = "MA Plot", ylim = c(-5, 5))
dev.off()

# 火山图
cat("生成火山图...\n")
res_df <- as.data.frame(res)
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up"
res_df$significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down"
res_df$significance <- factor(res_df$significance, levels = c("Down", "Not Significant", "Up"))

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw() +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "right") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

# 添加标签（top 10 基因）
top_genes <- head(res_df[order(res_df$padj), ], 10)
p_volcano <- p_volcano + geom_text_repel(
  data = top_genes,
  aes(label = rownames(top_genes)),
  size = 3,
  max.overlaps = 20
)

ggsave(file.path(output_dir, "volcano_plot.png"), p_volcano, width = 10, height = 8)

# PCA 图
cat("生成 PCA 图...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA Plot",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

ggsave(file.path(output_dir, "PCA_plot.png"), p_pca, width = 8, height = 6)

# 热图
cat("生成热图...\n")
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[top_var_genes, ]
mat <- mat - rowMeans(mat)

png(file.path(output_dir, "heatmap_top50.png"), width = 800, height = 1000)
pheatmap(
  mat,
  annotation_col = as.data.frame(colData(vsd)["condition"]),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "Top 50 Variable Genes"
)
dev.off()

# 样本距离热图
cat("生成样本距离热图...\n")
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

png(file.path(output_dir, "sample_distance_heatmap.png"), width = 800, height = 800)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  main = "Sample Distance Heatmap"
)
dev.off()

# 保存标准化表达矩阵
cat("保存标准化表达矩阵...\n")
write.csv(as.data.frame(assay(vsd)), file = file.path(output_dir, "normalized_counts.csv"))

# 保存 DESeq2 对象
save(dds, res, vsd, file = file.path(output_dir, "deseq2_results.RData"))

cat("==========================================\n")
cat("DESeq2 分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
