#!/usr/bin/env Rscript
# 差异 Peak 分析 (ATAC-seq)
# 用法: Rscript deseq2_peaks.R <peak_dir> <bam_dir> <sample_info> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
peak_dir <- ifelse(length(args) >= 1, args[1], "output/04_peak_calling")
bam_dir <- ifelse(length(args) >= 2, args[2], "output/03_alignment")
sample_info_file <- ifelse(length(args >= 3), args[3], "data/sample_info.txt")
output_dir <- ifelse(length(args) >= 4, args[4], "output/05_diff_peaks")

cat("==========================================\n")
cat("ATAC-seq 差异 Peak 分析\n")
cat("==========================================\n")
cat("Peak 目录:", peak_dir, "\n")
cat("BAM 目录:", bam_dir, "\n")
cat("样本信息:", sample_info_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==================== 准备 peak 区域 ====================

cat("\n读取合并的 peak 区域...\n")

# 读取合并的 peaks
merged_peaks <- read.table(file.path(dirname(peak_dir), "merged_peaks.bed"),
                           header = FALSE, stringsAsFactors = FALSE)
colnames(merged_peaks) <- c("chr", "start", "end")

# 创建 GRanges 对象
peaks_gr <- GRanges(
  seqnames = merged_peaks$chr,
  ranges = IRanges(start = merged_peaks$start + 1, end = merged_peaks$end)
)

cat("合并的 peak 数量:", length(peaks_gr), "\n")

# ==================== 计算 peak counts ====================

cat("\n计算每个样本的 peak counts...\n")

# 读取样本信息
sample_info <- read.table(sample_info_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# 获取 BAM 文件列表
bam_files <- list.files(bam_dir, pattern = "\\.final\\.bam$", recursive = TRUE, full.names = TRUE)
sample_names <- sapply(strsplit(basename(bam_files), "\\."), function(x) x[1])

cat("找到样本:", paste(sample_names, collapse = ", "), "\n")

# 使用 GenomicAlignments 计算 counts
library(GenomicAlignments)

count_matrix <- matrix(0, nrow = length(peaks_gr), ncol = length(bam_files))
colnames(count_matrix) <- sample_names
rownames(count_matrix) <- paste0("peak_", seq_len(length(peaks_gr)))

for (i in seq_along(bam_files)) {
  cat("  处理:", sample_names[i], "\n")

  # 读取 BAM 文件
  bam <- readGAlignments(bam_files[i], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE)))

  # 计算 overlaps
  overlaps <- findOverlaps(peaks_gr, bam, ignore.strand = TRUE)
  counts <- table(queryHits(overlaps))

  count_matrix[names(counts), i] <- as.numeric(counts)
}

# 保存 count matrix
write.csv(count_matrix, file = file.path(output_dir, "peak_counts.csv"))

# ==================== DESeq2 分析 ====================

cat("\n运行 DESeq2 分析...\n")

# 过滤低计数 peaks
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix_filtered <- count_matrix[keep, ]
cat("过滤后 peak 数量:", nrow(count_matrix_filtered), "\n")

# 确保样本顺序一致
common_samples <- intersect(colnames(count_matrix_filtered), rownames(sample_info))
count_matrix_filtered <- count_matrix_filtered[, common_samples]
sample_info <- sample_info[common_samples, , drop = FALSE]

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = sample_info,
  design = ~ condition
)

# 运行 DESeq2
dds <- DESeq(dds)

# 获取结果
res <- results(dds, alpha = 0.05)

# lfc shrinkage
res <- lfcShrink(dds, coef = 2, type = "apeglm")

# 排序
res_ordered <- res[order(res$padj), ]

# 添加基因组坐标
peak_coords <- data.frame(
  peak_id = rownames(res_ordered),
  chr = as.character(seqnames(peaks_gr[as.numeric(gsub("peak_", "", rownames(res_ordered)))])),
  start = start(peaks_gr[as.numeric(gsub("peak_", "", rownames(res_ordered)))]),
  end = end(peaks_gr[as.numeric(gsub("peak_", "", rownames(res_ordered)))])
)

# 合并结果
res_df <- cbind(peak_coords, as.data.frame(res_ordered))

# 保存所有结果
write.csv(res_df, file = file.path(output_dir, "diff_peaks_all.csv"), row.names = FALSE)

# 筛选显著差异 peaks
sig_peaks <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_peaks, file = file.path(output_dir, "diff_peaks_significant.csv"), row.names = FALSE)

cat("\n差异 Peak 统计:\n")
cat("  总 Peak 数:", nrow(res_df), "\n")
cat("  显著差异 Peak 数:", nrow(sig_peaks), "\n")
cat("  上调:", sum(sig_peaks$log2FoldChange > 0), "\n")
cat("  下调:", sum(sig_peaks$log2FoldChange < 0), "\n")

# ==================== 可视化 ====================

cat("\n生成可视化图表...\n")

# 火山图
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Open"
res_df$significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Closed"
res_df$significance <- factor(res_df$significance, levels = c("Closed", "Not Significant", "Open"))

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Closed" = "#2166AC", "Not Significant" = "grey70", "Open" = "#B2182B")) +
  theme_bw() +
  labs(
    title = "Differential Accessibility",
    subtitle = paste0("Open: ", sum(res_df$significance == "Open"), " | Closed: ", sum(res_df$significance == "Closed")),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

ggsave(file.path(output_dir, "volcano_plot.png"), p_volcano, width = 10, height = 8)

# MA 图
p_ma <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  scale_color_manual(values = c("Closed" = "#2166AC", "Not Significant" = "grey70", "Open" = "#B2182B")) +
  theme_bw() +
  labs(
    title = "MA Plot",
    x = "Mean of Normalized Counts",
    y = "log2 Fold Change"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "MA_plot.png"), p_ma, width = 10, height = 8)

# PCA 图
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
top_var_peaks <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[top_var_peaks, ]
mat <- mat - rowMeans(mat)

png(file.path(output_dir, "heatmap_top50.png"), width = 800, height = 1000)
pheatmap(
  mat,
  annotation_col = as.data.frame(colData(vsd)["condition"]),
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "Top 50 Variable Peaks"
)
dev.off()

# ==================== 导出显著 peaks 为 BED ====================

cat("\n导出显著 peaks...\n")

# Open peaks
open_peaks <- sig_peaks[sig_peaks$log2FoldChange > 0, c("chr", "start", "end", "peak_id", "log2FoldChange", "padj")]
open_peaks_gr <- GRanges(
  seqnames = open_peaks$chr,
  ranges = IRanges(start = open_peaks$start, end = open_peaks$end),
  name = open_peaks$peak_id,
  score = open_peaks$log2FoldChange
)
export.bed(open_peaks_gr, file.path(output_dir, "open_peaks.bed"))

# Closed peaks
closed_peaks <- sig_peaks[sig_peaks$log2FoldChange < 0, c("chr", "start", "end", "peak_id", "log2FoldChange", "padj")]
closed_peaks_gr <- GRanges(
  seqnames = closed_peaks$chr,
  ranges = IRanges(start = closed_peaks$start, end = closed_peaks$end),
  name = closed_peaks$peak_id,
  score = abs(closed_peaks$log2FoldChange)
)
export.bed(closed_peaks_gr, file.path(output_dir, "closed_peaks.bed"))

# 保存 DESeq2 对象
save(dds, res, vsd, peaks_gr, file = file.path(output_dir, "deseq2_peaks.RData"))

cat("\n==========================================\n")
cat("差异 Peak 分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
