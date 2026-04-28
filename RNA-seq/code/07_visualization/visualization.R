#!/usr/bin/env Rscript
# 数据可视化
# 用法: Rscript visualization.R <deseq2_results> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(VennDiagram)
  library(ggpubr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
deseq2_file <- ifelse(length(args) >= 1, args[1], "output/05_de_analysis/deseq2_results.RData")
output_dir <- ifelse(length(args) >= 2, args[2], "output/07_visualization")

cat("==========================================\n")
cat("数据可视化\n")
cat("==========================================\n")
cat("DESeq2 结果:", deseq2_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 加载 DESeq2 结果
cat("加载 DESeq2 结果...\n")
load(deseq2_file)

# 获取数据
normalized_counts <- as.data.frame(assay(vsd))
res_df <- as.data.frame(res)

# ==================== 差异表达基因可视化 ====================

cat("\n生成差异表达分析图表...\n")

# 火山图（增强版）
deg_df <- res_df %>%
  mutate(
    significance = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "Up",
      log2FoldChange < -1 & padj < 0.05 ~ "Down",
      TRUE ~ "Not Significant"
    )
  ) %>%
  mutate(significance = factor(significance, levels = c("Down", "Not Significant", "Up")))

p_volcano <- ggplot(deg_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "#2166AC", "Not Significant" = "#B2182B", "Up" = "#B2182B")) +
  theme_classic() +
  labs(
    title = "Volcano Plot",
    subtitle = paste0("Up: ", sum(deg_df$significance == "Up"), " | Down: ", sum(deg_df$significance == "Down")),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

# 标注 top 基因
top_genes <- deg_df %>%
  arrange(padj) %>%
  head(15)

p_volcano <- p_volcano +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.5
  )

ggsave(file.path(output_dir, "volcano_plot_enhanced.png"), p_volcano, width = 12, height = 10)

# ==================== 热图 ====================

cat("生成热图...\n")

# 筛选差异基因
sig_genes <- rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1))

if (length(sig_genes) > 0) {
  # 限制基因数量
  if (length(sig_genes) > 100) {
    sig_genes <- head(sig_genes[order(abs(res_df[sig_genes, "log2FoldChange"]), decreasing = TRUE)], 100)
  }

  heatmap_data <- normalized_counts[sig_genes, ]

  # 使用 ComplexHeatmap 绘制热图
  ha <- HeatmapAnnotation(
    Condition = colData(vsd)$condition,
    col = list(Condition = setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), unique(colData(vsd)$condition)))
  )

  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  ht <- Heatmap(
    heatmap_data,
    name = "Z-score",
    top_annotation = ha,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "Expression\n(Z-score)",
      at = c(-2, 0, 2),
      labels = c("Low", "Medium", "High")
    )
  )

  png(file.path(output_dir, "heatmap_DEG.png"), width = 1000, height = 1200, res = 120)
  draw(ht)
  dev.off()
}

# ==================== 样本关系图 ====================

cat("生成样本关系图...\n")

# PCA 图（增强版）
pca_res <- prcomp(t(assay(vsd)), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Sample = rownames(pca_res$x),
  Condition = colData(vsd)$condition
)

var_explained <- summary(pca_res)$importance[2, ] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Condition)) +
  theme_classic() +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "PCA_enhanced.png"), p_pca, width = 10, height = 8)

# 样本聚类树
sample_dist <- dist(t(assay(vsd)))
hc <- hclust(sample_dist)

png(file.path(output_dir, "sample_clustering.png"), width = 800, height = 600)
plot(hc, main = "Sample Clustering Dendrogram", xlab = "Samples", sub = "")
dev.off()

# ==================== 箱线图 ====================

cat("生成表达分布图...\n")

# 表达分布箱线图
expr_long <- normalized_counts %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(data.frame(Sample = rownames(colData(vsd)), Condition = colData(vsd)$condition), by = "Sample")

p_boxplot <- ggplot(expr_long, aes(x = Sample, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(title = "Expression Distribution", y = "Normalized Expression") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "expression_boxplot.png"), p_boxplot, width = 12, height = 6)

# ==================== 小提琴图 ====================

cat("生成小提琴图...\n")

p_violin <- ggplot(expr_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(title = "Expression Distribution by Condition", y = "Normalized Expression") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(file.path(output_dir, "expression_violin.png"), p_violin, width = 10, height = 6)

# ==================== 特定基因表达图 ====================

cat("生成特定基因表达图...\n")

# 选择 top 10 差异基因
top_genes_plot <- head(rownames(res_df[order(res_df$padj), ]), 10)

for (gene in top_genes_plot) {
  gene_expr <- data.frame(
    Sample = colnames(normalized_counts),
    Expression = as.numeric(normalized_counts[gene, ]),
    Condition = colData(vsd)$condition
  )

  p <- ggplot(gene_expr, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    theme_classic() +
    labs(
      title = paste0("Expression of ", gene),
      y = "Normalized Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )

  ggsave(file.path(output_dir, paste0("gene_", gene, ".png")), p, width = 6, height = 5)
}

# ==================== 相关性热图 ====================

cat("生成相关性热图...\n")

cor_matrix <- cor(assay(vsd), method = "pearson")

p_corr <- pheatmap(
  cor_matrix,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Sample Correlation Matrix",
  fontsize_number = 8
)

png(file.path(output_dir, "correlation_heatmap.png"), width = 800, height = 800)
print(p_corr)
dev.off()

# ==================== MA 图 ====================

cat("生成 MA 图...\n")

ma_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  scale_color_manual(values = c("Down" = "#2166AC", "NS" = "grey70", "Up" = "#B2182B")) +
  theme_classic() +
  labs(
    title = "MA Plot",
    x = "Mean of Normalized Counts",
    y = "log2 Fold Change"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "MA_plot_enhanced.png"), p_ma, width = 10, height = 8)

# ==================== 生成汇总报告 ====================

cat("\n生成分析汇总...\n")

summary_df <- data.frame(
  Metric = c(
    "Total Genes Analyzed",
    "Up-regulated Genes (log2FC > 1, padj < 0.05)",
    "Down-regulated Genes (log2FC < -1, padj < 0.05)",
    "Total Samples",
    "Number of Conditions"
  ),
  Value = c(
    nrow(res_df),
    sum(ma_df$significance == "Up"),
    sum(ma_df$significance == "Down"),
    ncol(normalized_counts),
    length(unique(colData(vsd)$condition))
  )
)

write.csv(summary_df, file = file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

cat("==========================================\n")
cat("可视化分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
