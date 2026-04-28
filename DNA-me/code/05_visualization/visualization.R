#!/usr/bin/env Rscript
# DNA 甲基化可视化
# 用法: Rscript visualization.R <dmr_results> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(pheatmap)
  library(methylKit)
  library(bsseq)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
dmr_file <- ifelse(length(args) >= 1, args[1], "output/03_dmr_analysis/dmr_analysis_results.RData")
output_dir <- ifelse(length(args) >= 2, args[2], "output/05_visualization")

cat("==========================================\n")
cat("DNA 甲基化可视化\n")
cat("==========================================\n")
cat("DMR 结果:", dmr_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 加载数据
cat("\n加载数据...\n")
load(dmr_file)
load("output/02_calculation/methylation_objects.RData")

# 读取 DMR 和 DML 数据
dmr_results <- read.csv("output/03_dmr_analysis/DMR_all.csv")
dml_results <- read.csv("output/03_dmr_analysis/DML_all.csv")
meth_matrix <- read.csv("output/02_calculation/methylation_matrix.csv")

# ==================== 甲基化水平分布 ====================

cat("\n生成甲基化水平分布图...\n")

# 密度图
meth_long <- meth_matrix %>%
  select(-chr, -start, -end, -strand) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Methylation") %>%
  filter(!is.na(Methylation))

p_density <- ggplot(meth_long, aes(x = Methylation, fill = Sample, color = Sample)) +
  geom_density(alpha = 0.3) +
  theme_classic() +
  labs(
    title = "Methylation Level Distribution",
    x = "Methylation Level (%)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "methylation_distribution.png"), p_density, width = 12, height = 6)

# 箱线图
sample_info <- read.table("data/sample_info.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

p_boxplot <- meth_long %>%
  left_join(sample_info, by = c("Sample" = "sample_id")) %>%
  ggplot(aes(x = Sample, y = Methylation, fill = condition)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "Methylation Level by Sample",
    x = "Sample",
    y = "Methylation Level (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "methylation_boxplot.png"), p_boxplot, width = 12, height = 6)

# ==================== DML 火山图 ====================

cat("\n生成 DML 火山图...\n")

dml_results$significance <- "Not Significant"
dml_results$significance[dml_results$diff > 0.25 & dml_results$fdr < 0.05] <- "Hyper"
dml_results$significance[dml_results$diff < -0.25 & dml_results$fdr < 0.05] <- "Hypo"
dml_results$significance <- factor(dml_results$significance, levels = c("Hypo", "Not Significant", "Hyper"))

p_volcano <- ggplot(dml_results, aes(x = diff * 100, y = -log10(fdr), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Hypo" = "#2166AC", "Not Significant" = "grey70", "Hyper" = "#B2182B")) +
  theme_classic() +
  labs(
    title = "Differential Methylation Loci",
    subtitle = paste0("Hyper: ", sum(dml_results$significance == "Hyper"),
                      " | Hypo: ", sum(dml_results$significance == "Hypo")),
    x = "Methylation Difference (%)",
    y = "-log10(FDR)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

ggsave(file.path(output_dir, "DML_volcano_enhanced.png"), p_volcano, width = 12, height = 10)

# ==================== DMR 热图 ====================

cat("\n生成 DMR 热图...\n")

# 获取 DMR 区域的甲基化值
dmr_gr <- GRanges(
  seqnames = dmr_results$chr,
  ranges = IRanges(start = dmr_results$start, end = dmr_results$end)
)

# 筛选显著的 DMR
sig_dmrs <- dmr_results[abs(dmr_results$mean_diff) > 0.2 & dmr_results$nCG >= 5, ]

if (nrow(sig_dmrs) > 0) {
  # 取 top DMR
  top_dmrs <- head(sig_dmrs[order(abs(sig_dmrs$areaStat), decreasing = TRUE), ], 50)

  # 创建甲基化矩阵
  dmr_meth <- data.frame(
    DMR = paste0("DMR_", seq_len(nrow(top_dmrs))),
    Direction = top_dmrs$direction
  )

  # 添加每个样本的平均甲基化水平
  for (sample in colnames(meth_matrix)[5:ncol(meth_matrix)]) {
    # 计算每个 DMR 区域的平均甲基化
    dmr_meth[[sample]] <- sapply(seq_len(nrow(top_dmrs)), function(i) {
      dmr <- top_dmrs[i, ]
      sites <- meth_matrix %>%
        filter(chr == dmr$chr, start >= dmr$start, end <= dmr$end)
      if (nrow(sites) > 0) {
        mean(sites[[sample]], na.rm = TRUE)
      } else {
        NA
      }
    })
  }

  # 绘制热图
  mat <- as.matrix(dmr_meth[, 3:ncol(dmr_meth)])
  rownames(mat) <- dmr_meth$DMR

  # 去除 NA
  mat[is.na(mat)] <- 50

  # 缩放
  mat_scaled <- t(scale(t(mat)))

  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  ha <- rowAnnotation(
    Direction = dmr_meth$Direction,
    col = list(Direction = c("Hyper" = "#B2182B", "Hypo" = "#2166AC"))
  )

  ht <- Heatmap(
    mat_scaled,
    name = "Z-score",
    left_annotation = ha,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = col_fun,
    column_names_gp = gpar(fontsize = 10),
    row_title = "DMRs",
    column_title = "Samples"
  )

  png(file.path(output_dir, "DMR_heatmap.png"), width = 1000, height = 1200, res = 120)
  draw(ht)
  dev.off()
}

# ==================== DMR 基因组分布 ====================

cat("\n生成 DMR 基因组分布图...\n")

# 染色体分布
chr_counts <- dmr_results %>%
  group_by(chr) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

p_chr <- ggplot(head(chr_counts, 20), aes(x = reorder(chr, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#3182BD") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "DMR Distribution Across Chromosomes",
    x = "Chromosome",
    y = "Number of DMRs"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "DMR_chromosome_distribution.png"), p_chr, width = 10, height = 8)

# DMR 长度分布
p_length <- ggplot(dmr_results, aes(x = length)) +
  geom_histogram(bins = 50, fill = "#3182BD", color = "white") +
  theme_classic() +
  labs(
    title = "DMR Length Distribution",
    x = "DMR Length (bp)",
    y = "Count"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "DMR_length_distribution.png"), p_length, width = 10, height = 6)

# 方向分布
p_direction <- ggplot(dmr_results, aes(x = direction, fill = direction)) +
  geom_bar() +
  scale_fill_manual(values = c("Hyper" = "#B2182B", "Hypo" = "#2166AC")) +
  theme_classic() +
  labs(
    title = "DMR Direction Distribution",
    x = "Direction",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "DMR_direction_distribution.png"), p_direction, width = 6, height = 6)

# ==================== 样本相关性 ====================

cat("\n生成样本相关性分析...\n")

# 相关性矩阵
cor_matrix <- cor(as.matrix(meth_matrix[, 5:ncol(meth_matrix)]), use = "complete.obs", method = "pearson")

p_corr <- pheatmap(
  cor_matrix,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Sample Correlation Matrix",
  fontsize_number = 8
)

png(file.path(output_dir, "sample_correlation_heatmap.png"), width = 800, height = 800)
print(p_corr)
dev.off()

# ==================== PCA 分析 ====================

cat("\n生成 PCA 分析...\n")

# 对甲基化矩阵进行 PCA
meth_mat <- as.matrix(meth_matrix[, 5:ncol(meth_matrix)])
meth_mat[is.na(meth_mat)] <- 0

pca_res <- prcomp(t(meth_mat), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Sample = colnames(meth_mat)
) %>%
  left_join(sample_info, by = c("Sample" = "sample_id"))

var_explained <- summary(pca_res)$importance[2, ] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(file.path(output_dir, "PCA_plot.png"), p_pca, width = 10, height = 8)

# ==================== CpG 密度与甲基化关系 ====================

cat("\n生成 CpG 密度分析...\n")

cpg_density <- read.csv("output/02_calculation/cpg_density_1kb.csv")

if (nrow(cpg_density) > 0) {
  # CpG 密度分布
  p_cpg <- ggplot(cpg_density, aes(x = .)) +
    geom_histogram(bins = 50, fill = "#3182BD") +
    theme_classic() +
    labs(
      title = "CpG Density Distribution (1kb windows)",
      x = "CpG Count",
      y = "Frequency"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, "cpg_density_distribution.png"), p_cpg, width = 10, height = 6)
}

# ==================== 汇总统计 ====================

cat("\n生成汇总统计...\n")

summary_stats <- data.frame(
  Metric = c(
    "Total DMRs",
    "Hyper-methylated DMRs",
    "Hypo-methylated DMRs",
    "Mean DMR Length (bp)",
    "Total DMLs",
    "Significant DMLs (FDR < 0.05, |diff| > 25%)"
  ),
  Value = c(
    nrow(dmr_results),
    sum(dmr_results$direction == "Hyper"),
    sum(dmr_results$direction == "Hypo"),
    round(mean(dmr_results$length), 1),
    nrow(dml_results),
    sum(dml_results$fdr < 0.05 & abs(dml_results$diff) > 0.25)
  )
)

write.csv(summary_stats, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

cat("\n==========================================\n")
cat("可视化分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
