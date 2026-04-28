#!/usr/bin/env Rscript
# DSS 差异甲基化区域 (DMR) 分析
# 用法: Rscript dss_dmr_analysis.R <methylation_file> <sample_info> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(annotatr)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
methylation_file <- ifelse(length(args) >= 1, args[1], "output/02_calculation/methylation_objects.RData")
sample_info_file <- ifelse(length(args) >= 2, args[2], "data/sample_info.txt")
output_dir <- ifelse(length(args) >= 3, args[3], "output/03_dmr_analysis")

cat("==========================================\n")
cat("DSS 差异甲基化区域分析\n")
cat("==========================================\n")
cat("甲基化数据:", methylation_file, "\n")
cat("样本信息:", sample_info_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==================== 加载数据 ====================

cat("\n加载甲基化数据...\n")
load(methylation_file)

# 读取样本信息
sample_info <- read.table(sample_info_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# 获取处理组和对照组
treatment_samples <- sample_info$sample_id[sample_info$condition == "treatment"]
control_samples <- sample_info$sample_id[sample_info$condition == "control"]

cat("处理组样本:", paste(treatment_samples, collapse = ", "), "\n")
cat("对照组样本:", paste(control_samples, collapse = ", "), "\n")

# ==================== 准备 DSS 输入 ====================

cat("\n准备 DSS 分析数据...\n")

# 从 methylKit 对象提取数据
meth_matrix <- read.csv("output/02_calculation/methylation_matrix.csv")

# 创建 BSseq 对象
# 需要甲基化计数和总计数矩阵
createDSSInput <- function(filtered_meth, sample_names) {
  # 初始化矩阵
  M_list <- list()
  Cov_list <- list()

  for (i in seq_along(filtered_meth)) {
    obj <- filtered_meth[[i]]
    data <- getData(obj)

    # 甲基化 reads 数
    M_list[[i]] <- data$numCs
    # 总覆盖度
    Cov_list[[i]] <- data$coverage
  }

  # 合并为矩阵
  M <- do.call(cbind, M_list)
  Cov <- do.call(cbind, Cov_list)
  colnames(M) <- sample_names
  colnames(Cov) <- sample_names

  # 创建 GRanges
  gr <- GRanges(
    seqnames = getData(filtered_meth[[1]])$chr,
    ranges = IRanges(start = getData(filtered_meth[[1]])$start,
                     end = getData(filtered_meth[[1]])$end)
  )

  return(list(M = M, Cov = Cov, gr = gr))
}

dss_input <- createDSSInput(filtered_meth, names(filtered_meth))

# 创建 BSseq 对象
BSobj <- BSseq(
  gr = dss_input$gr,
  M = dss_input$M,
  Cov = dss_input$Cov,
  sampleNames = names(filtered_meth)
)

cat("BSseq 对象创建完成\n")
cat("位点数:", length(BSobj), "\n")

# ==================== DML 检测 ====================

cat("\n检测差异甲基化位点 (DML)...\n")

# 设置设计矩阵
design <- model.matrix(~ sample_info$condition)
colnames(design) <- c("Intercept", "Treatment")

# DML 检测
dml_fit <- DMLfit.multiFactor(BSobj, design = design)
dml_test <- DMLtest.multiFactor(dml_fit, coef = "Treatment")

# 保存 DML 结果
dml_results <- data.frame(
  chr = as.character(seqnames(dml_test)),
  start = start(dml_test),
  end = end(dml_test),
  mu1 = dml_test$mu1,
  mu2 = dml_test$mu2,
  diff = dml_test$diff,
  stat = dml_test$stat,
  pvalue = dml_test$pvals,
  fdr = p.adjust(dml_test$pvals, method = "BH")
)

# 排序
dml_results <- dml_results[order(dml_results$fdr), ]

# 保存所有 DML
write.csv(dml_results, file = file.path(output_dir, "DML_all.csv"), row.names = FALSE)

# 筛选显著 DML
sig_dml <- subset(dml_results, fdr < 0.05 & abs(diff) > 0.25)
write.csv(sig_dml, file = file.path(output_dir, "DML_significant.csv"), row.names = FALSE)

cat("\nDML 统计:\n")
cat("  总位点数:", nrow(dml_results), "\n")
cat("  显著 DML 数:", nrow(sig_dml), "\n")
cat("  高甲基化:", sum(sig_dml$diff > 0), "\n")
cat("  低甲基化:", sum(sig_dml$diff < 0), "\n")

# ==================== DMR 检测 ====================

cat("\n检测差异甲基化区域 (DMR)...\n")

# 使用 DMR 检测
dmrs <- callDMR(dml_test,
                p.threshold = 0.001,
                minlen = 100,
                minCG = 3,
                dis.merge = 100,
                pct.sig = 0.5)

# 转换为数据框
dmr_results <- data.frame(
  chr = dmrs$chr,
  start = dmrs$start,
  end = dmrs$end,
  length = dmrs$length,
  nCG = dmrs$nCG,
  mean_diff = dmrs$mean_diff,
  areaStat = dmrs$areaStat,
  direction = ifelse(dmrs$mean_diff > 0, "Hyper", "Hypo")
)

# 按区域统计排序
dmr_results <- dmr_results[order(abs(dmr_results$areaStat), decreasing = TRUE), ]

# 保存 DMR 结果
write.csv(dmr_results, file = file.path(output_dir, "DMR_all.csv"), row.names = FALSE)

# 筛选显著 DMR
sig_dmrs <- subset(dmr_results, abs(mean_diff) > 0.2 & nCG >= 5)
write.csv(sig_dmrs, file = file.path(output_dir, "DMR_significant.csv"), row.names = FALSE)

cat("\nDMR 统计:\n")
cat("  总 DMR 数:", nrow(dmr_results), "\n")
cat("  显著 DMR 数:", nrow(sig_dmrs), "\n")
cat("  高甲基化区域:", sum(sig_dmrs$direction == "Hyper"), "\n")
cat("  低甲基化区域:", sum(sig_dmrs$direction == "Hypo"), "\n")

# ==================== 导出 BED 文件 ====================

cat("\n导出 BED 文件...\n")

# DMR BED
dmr_gr <- GRanges(
  seqnames = dmr_results$chr,
  ranges = IRanges(start = dmr_results$start, end = dmr_results$end),
  name = paste0("DMR_", seq_len(nrow(dmr_results))),
  score = abs(dmr_results$mean_diff)
)

export.bed(dmr_gr, file.path(output_dir, "DMR_all.bed"))

# 高甲基化区域
hyper_gr <- dmr_gr[dmr_results$direction == "Hyper"]
export.bed(hyper_gr, file.path(output_dir, "DMR_hyper.bed"))

# 低甲基化区域
hypo_gr <- dmr_gr[dmr_results$direction == "Hypo"]
export.bed(hypo_gr, file.path(output_dir, "DMR_hypo.bed"))

# DML BED
dml_gr <- GRanges(
  seqnames = sig_dml$chr,
  ranges = IRanges(start = sig_dml$start, end = sig_dml$end),
  name = paste0("DML_", seq_len(nrow(sig_dml))),
  score = abs(sig_dml$diff)
)

export.bed(dml_gr, file.path(output_dir, "DML_significant.bed"))

# ==================== 可视化 ====================

cat("\n生成可视化图表...\n")

# 火山图
p_volcano <- ggplot(dml_results, aes(x = diff * 100, y = -log10(fdr))) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Hypo" = "#2166AC", "NotSig" = "grey70", "Hyper" = "#B2182B")) +
  theme_classic() +
  labs(
    title = "Differential Methylation Loci",
    x = "Methylation Difference (%)",
    y = "-log10(FDR)"
  ) +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "DML_volcano.png"), p_volcano, width = 10, height = 8)

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

# DMR 分布条形图
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

ggsave(file.path(output_dir, "DMR_direction.png"), p_direction, width = 6, height = 6)

# ==================== 保存对象 ====================

save(dml_fit, dml_test, dml_results, dmr_results, sig_dml, sig_dmrs,
     file = file.path(output_dir, "dmr_analysis_results.RData"))

cat("\n==========================================\n")
cat("DMR 分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
