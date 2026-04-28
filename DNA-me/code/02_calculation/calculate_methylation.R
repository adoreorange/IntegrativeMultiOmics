#!/usr/bin/env Rscript
# DNA 甲基化水平计算
# 用法: Rscript calculate_methylation.R <preprocess_dir> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(methylKit)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(data.table)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
preprocess_dir <- ifelse(length(args) >= 1, args[1], "output/01_preprocess")
output_dir <- ifelse(length(args) >= 2, args[2], "output/02_calculation")
sample_info_file <- ifelse(length(args) >= 3, args[3], "data/sample_info.txt")

cat("==========================================\n")
cat("DNA Methylation 计算\n")
cat("==========================================\n")
cat("预处理目录:", preprocess_dir, "\n")
cat("输出目录:", output_dir, "\n")
cat("样本信息:", sample_info_file, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==================== 读取样本信息 ====================

cat("\n读取样本信息...\n")
sample_info <- read.table(sample_info_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cat("样本数量:", nrow(sample_info), "\n")

# ==================== 读取甲基化数据 ====================

cat("\n读取甲基化数据...\n")

# 查找 CX_report 文件
cx_files <- list.files(preprocess_dir, pattern = "CX_report.txt$", recursive = TRUE, full.names = TRUE)

if (length(cx_files) == 0) {
  # 尝试其他格式
  cx_files <- list.files(preprocess_dir, pattern = "_bismark_bt2.*.deduplicated.CX_report.txt$", recursive = TRUE, full.names = TRUE)
}

if (length(cx_files) == 0) {
  stop("未找到 CX_report 文件，请检查预处理步骤")
}

cat("找到文件:", length(cx_files), "\n")

# 使用 methylKit 读取数据
sample_names <- sapply(basename(cx_files), function(x) {
  gsub("_bismark_bt2.*", "", gsub(".deduplicated.CX_report.txt", "", x))
})

# 读取 methylKit 对象
meth_objects <- list()

for (i in seq_along(cx_files)) {
  cat("  处理:", sample_names[i], "\n")

  # 获取样本信息
  sample_row <- sample_info[grepl(sample_names[i], sample_info$sample_id), ]

  if (nrow(sample_row) == 0) {
    warning(paste0("样本 ", sample_names[i], " 未在样本信息中找到，使用默认设置"))
    condition <- "control"
  } else {
    condition <- as.character(sample_row$condition[1])
  }

  # 读取数据
  meth_objects[[i]] <- methRead(
    cx_files[i],
    sample.id = sample_names[i],
    assembly = "hg38",
    pipeline = "bismarkCytosineReport",
    context = "CpG",
    treatment = ifelse(condition == "treatment", 1, 0)
  )
}

names(meth_objects) <- sample_names

# ==================== 过滤低覆盖度位点 ====================

cat("\n过滤低覆盖度位点...\n")

filtered_meth <- list()
for (i in seq_along(meth_objects)) {
  cat("  过滤:", sample_names[i], "\n")

  # 过滤覆盖度 < 10 或 > 99.9 百分位的位点
  filtered_meth[[i]] <- filterCoverage(
    meth_objects[[i]],
    lo.count = 10,
    lo.perc = NULL,
    hi.perc = 99.9
  )
}
names(filtered_meth) <- sample_names

# ==================== 合并样本 ====================

cat("\n合并样本数据...\n")

# 找到所有样本共有的 CpG 位点
merged_meth <- unite(filtered_meth, destrand = FALSE)

cat("共有 CpG 位点数:", nrow(merged_meth), "\n")

# 保存合并数据
save(merged_meth, file = file.path(output_dir, "merged_methylation.RData"))

# ==================== 计算甲基化水平统计 ====================

cat("\n计算甲基化统计...\n")

# 获取甲基化百分比矩阵
meth_matrix <- getData(merged_meth)

# 每个样本的统计
sample_stats <- data.frame(
  Sample = sample_names,
  Total_CpGs = ncol(meth_matrix) - 4,  # 减去 chr, start, end, strand 列
  Mean_Methylation = colMeans(meth_matrix[, grep("freqC", colnames(meth_matrix))] /
                                (meth_matrix[, grep("freqC", colnames(meth_matrix))] +
                                   meth_matrix[, grep("freqT", colnames(meth_matrix))]), na.rm = TRUE),
  Median_Coverage = apply(meth_matrix[, grep("coverage", colnames(meth_matrix))], 2, median, na.rm = TRUE)
)

write.csv(sample_stats, file = file.path(output_dir, "sample_methylation_stats.csv"), row.names = FALSE)

# ==================== 生成基因组范围的甲基化矩阵 ====================

cat("\n生成甲基化矩阵...\n")

# 提取甲基化百分比
meth_pct <- getData(merged_meth)

# 创建 GRanges 对象
meth_gr <- GRanges(
  seqnames = meth_pct$chr,
  ranges = IRanges(start = meth_pct$start, end = meth_pct$end),
  strand = meth_pct$strand
)

# 提取每个样本的甲基化百分比
methylation_matrix <- data.frame(
  chr = meth_pct$chr,
  start = meth_pct$start,
  end = meth_pct$end,
  strand = meth_pct$strand
)

for (i in seq_along(filtered_meth)) {
  sample_name <- sample_names[i]
  freqC_col <- paste0(sample_name, ".freqC")
  coverage_col <- paste0(sample_name, ".coverage")

  # 计算甲基化百分比
  if (freqC_col %in% colnames(meth_pct)) {
    methylation_matrix[[sample_name]] <- meth_pct[[freqC_col]] / meth_pct[[coverage_col]] * 100
  }
}

# 保存甲基化矩阵
write.csv(methylation_matrix, file = file.path(output_dir, "methylation_matrix.csv"), row.names = FALSE)

# ==================== 计算 CpG 密度 ====================

cat("\n计算 CpG 密度...\n")

# 使用 tileMethylCounts 计算 CpG 密度
tiled_meth <- tileMethylCounts(merged_meth, win.size = 1000, step.size = 1000)

# 获取 CpG 密度
cpg_density <- getData(tiled_meth)

# 保存
write.csv(cpg_density, file = file.path(output_dir, "cpg_density_1kb.csv"), row.names = FALSE)

# ==================== 染色体分布统计 ====================

cat("\n计算染色体分布...\n")

chr_dist <- meth_pct %>%
  group_by(chr) %>%
  summarise(
    CpG_count = n(),
    Mean_methylation = mean(get(grep("freqC", colnames(.), value = TRUE)[1]) /
                              (get(grep("freqC", colnames(.), value = TRUE)[1]) +
                                 get(grep("freqT", colnames(.), value = TRUE)[1])) * 100, na.rm = TRUE)
  ) %>%
  arrange(desc(CpG_count))

write.csv(chr_dist, file = file.path(output_dir, "chromosome_distribution.csv"), row.names = FALSE)

# ==================== 保存对象 ====================

save(filtered_meth, meth_objects, sample_stats, methylation_matrix,
     file = file.path(output_dir, "methylation_objects.RData"))

cat("\n==========================================\n")
cat("甲基化计算完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
