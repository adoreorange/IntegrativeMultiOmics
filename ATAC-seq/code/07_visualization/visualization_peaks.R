#!/usr/bin/env Rscript
# ATAC-seq 可视化
# 用法: Rscript visualization_peaks.R <deseq2_results> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
deseq2_file <- ifelse(length(args) >= 1, args[1], "output/05_diff_peaks/deseq2_peaks.RData")
output_dir <- ifelse(length(args) >= 2, args[2], "output/07_visualization")

cat("==========================================\n")
cat("ATAC-seq 可视化分析\n")
cat("==========================================\n")
cat("DESeq2 结果:", deseq2_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 加载 DESeq2 结果
cat("\n加载 DESeq2 结果...\n")
load(deseq2_file)

# 获取数据
normalized_counts <- as.data.frame(assay(vsd))
res_df <- read.csv("output/05_diff_peaks/diff_peaks_significant.csv", stringsAsFactors = FALSE)

# ==================== 染色质开放性热图 ====================

cat("\n生成染色质开放性热图...\n")

# 准备数据
sig_peaks <- res_df[order(res_df$padj), ]
if (nrow(sig_peaks) > 100) {
  sig_peaks <- head(sig_peaks, 100)
}

# 获取 counts
peak_ids <- sig_peaks$peak_id
heatmap_data <- normalized_counts[peak_ids, , drop = FALSE]

# 复杂热图
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

ha <- HeatmapAnnotation(
  Condition = colData(vsd)$condition,
  col = list(Condition = setNames(c("#E41A1C", "#377EB8", "#4DAF4A"), unique(colData(vsd)$condition)))
)

ht <- Heatmap(
  as.matrix(heatmap_data),
  name = "Z-score",
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

png(file.path(output_dir, "chromatin_accessibility_heatmap.png"), width = 1000, height = 1200, res = 120)
draw(ht)
dev.off()

# ==================== 火山图 (增强版) ====================

cat("生成增强版火山图...\n")

res_all <- read.csv("output/05_diff_peaks/diff_peaks_all.csv", stringsAsFactors = FALSE)

res_all$significance <- "Not Significant"
res_all$significance[res_all$log2FoldChange > 1 & res_all$padj < 0.05] <- "Open"
res_all$significance[res_all$log2FoldChange < -1 & res_all$padj < 0.05] <- "Closed"
res_all$significance <- factor(res_all$significance, levels = c("Closed", "Not Significant", "Open"))

p_volcano <- ggplot(res_all, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Closed" = "#2166AC", "Not Significant" = "grey70", "Open" = "#B2182B")) +
  theme_classic() +
  labs(
    title = "Differential Chromatin Accessibility",
    subtitle = paste0("Open: ", sum(res_all$significance == "Open"),
                      " | Closed: ", sum(res_all$significance == "Closed")),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

ggsave(file.path(output_dir, "volcano_plot_enhanced.png"), p_volcano, width = 12, height = 10)

# ==================== Peak 分布图 ====================

cat("生成 Peak 基因组分布图...\n")

# 创建 GRanges
peaks_gr <- GRanges(
  seqnames = res_df$chr,
  ranges = IRanges(start = res_df$start, end = res_df$end)
)

# 染色体分布
chr_counts <- as.data.frame(table(seqnames(peaks_gr)))
colnames(chr_counts) <- c("Chromosome", "Count")
chr_counts <- chr_counts[order(chr_counts$Count, decreasing = TRUE), ]

p_chr <- ggplot(head(chr_counts, 20), aes(x = reorder(Chromosome, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#3182BD") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Peak Distribution Across Chromosomes",
    x = "Chromosome",
    y = "Number of Peaks"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "peak_chromosome_distribution.png"), p_chr, width = 10, height = 8)

# Peak 长度分布
peak_lengths <- width(peaks_gr)

p_length <- ggplot(data.frame(Length = peak_lengths), aes(x = Length)) +
  geom_histogram(bins = 50, fill = "#3182BD", color = "white") +
  theme_classic() +
  labs(
    title = "Peak Length Distribution",
    x = "Peak Length (bp)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  scale_x_log10()

ggsave(file.path(output_dir, "peak_length_distribution.png"), p_length, width = 10, height = 6)

# ==================== TSS 热图 ====================

cat("生成 TSS 区域热图...\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 获取 TSS
promoters <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

# 计算 peak 与 TSS 的关系
if (file.exists("output/04_peak_calling/merged_peaks.bed")) {
  all_peaks <- read.table("output/04_peak_calling/merged_peaks.bed", stringsAsFactors = FALSE)
  colnames(all_peaks) <- c("chr", "start", "end")

  all_peaks_gr <- GRanges(
    seqnames = all_peaks$chr,
    ranges = IRanges(start = all_peaks$start + 1, end = all_peaks$end)
  )

  # Peak 注释
  peak_anno <- annotatePeak(all_peaks_gr, tssRegion = c(-3000, 3000),
                            TxDb = txdb, annoDb = "org.Hs.eg.db")

  # 保存
  write.csv(as.data.frame(peak_anno), file.path(output_dir, "peak_TSS_annotation.csv"))

  # TSS 分布图
  p_tss <- plotDistToTSS(peak_anno) +
    ggtitle("Peak Distribution Relative to TSS") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, "peak_TSS_distribution.png"), p_tss, width = 10, height = 6)
}

# ==================== 样本相关性分析 ====================

cat("生成样本相关性分析...\n")

# 相关性矩阵
cor_matrix <- cor(as.matrix(normalized_counts), method = "pearson")

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

cat("生成 PCA 分析...\n")

pca_res <- prcomp(t(as.matrix(normalized_counts)), scale. = TRUE)
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

# ==================== FRiP 分析 ====================

cat("生成 FRiP 分析...\n")

if (file.exists("output/04_peak_calling/frip_scores.csv")) {
  frip <- read.csv("output/04_peak_calling/frip_scores.csv")

  p_frip <- ggplot(frip, aes(x = Sample, y = FRiP, fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(
      title = "FRiP Scores (Fraction of Reads in Peaks)",
      subtitle = "Dashed line indicates recommended minimum (0.3)",
      x = "Sample",
      y = "FRiP Score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, 1)

  ggsave(file.path(output_dir, "FRiP_scores.png"), p_frip, width = 10, height = 6)
}

# ==================== 汇总统计 ====================

cat("\n生成汇总统计...\n")

summary_stats <- data.frame(
  Metric = c(
    "Total Significant Peaks",
    "Open Peaks (log2FC > 1)",
    "Closed Peaks (log2FC < -1)",
    "Total Samples",
    "Number of Conditions"
  ),
  Value = c(
    nrow(res_df),
    sum(res_df$log2FoldChange > 0),
    sum(res_df$log2FoldChange < 0),
    ncol(normalized_counts),
    length(unique(colData(vsd)$condition))
  )
)

write.csv(summary_stats, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

cat("\n==========================================\n")
cat("可视化分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
