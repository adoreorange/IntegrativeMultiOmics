#!/usr/bin/env Rscript
# Peak 富集分析 (ATAC-seq)
# 用法: Rscript enrichment_peaks.R <diff_peaks_file> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(enrichplot)
  library(ggplot2)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
diff_peaks_file <- ifelse(length(args) >= 1, args[1], "output/05_diff_peaks/diff_peaks_significant.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "output/06_enrichment")

# 物种设置（可根据需要修改）
org_db <- org.Hs.eg.db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat("==========================================\n")
cat("ATAC-seq Peak 富集分析\n")
cat("==========================================\n")
cat("差异 Peak 文件:", diff_peaks_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取差异 peaks
cat("\n读取差异 peaks...\n")
diff_peaks <- read.csv(diff_peaks_file, stringsAsFactors = FALSE)

# 创建 GRanges 对象
open_peaks <- diff_peaks[diff_peaks$log2FoldChange > 0, ]
closed_peaks <- diff_peaks[diff_peaks$log2FoldChange < 0, ]

cat("Open peaks:", nrow(open_peaks), "\n")
cat("Closed peaks:", nrow(closed_peaks), "\n")

# 创建 GRanges
makeGRanges <- function(df) {
  if (nrow(df) == 0) return(NULL)
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end),
    name = df$peak_id
  )
}

open_gr <- makeGRanges(open_peaks)
closed_gr <- makeGRanges(closed_peaks)
all_gr <- makeGRanges(diff_peaks)

# ==================== Peak 注释 ====================

cat("\n进行 Peak 基因组注释...\n")

# 注释函数
annotatePeaks <- function(peaks, name) {
  if (is.null(peaks)) return(NULL)

  peak_anno <- annotatePeak(
    peaks,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )

  # 保存注释结果
  write.csv(as.data.frame(peak_anno), file.path(output_dir, paste0(name, "_annotation.csv")))

  # 绑定区域可视化
  p <- plotAnnoBar(peak_anno) +
    ggtitle(paste0(name, " - Genomic Annotation")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, paste0(name, "_annotation_bar.png")), p, width = 10, height = 6)

  # TSS 距离分布
  p_dist <- plotDistToTSS(peak_anno) +
    ggtitle(paste0(name, " - Distance to TSS")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, paste0(name, "_dist_to_tss.png")), p_dist, width = 10, height = 6)

  return(peak_anno)
}

# 注释所有 peaks
open_anno <- annotatePeaks(open_gr, "open_peaks")
closed_anno <- annotatePeaks(closed_gr, "closed_peaks")
all_anno <- annotatePeaks(all_gr, "all_peaks")

# ==================== 功能富集分析 ====================cat("\n进行功能富集分析...\n")

# 获取 peak 相关基因
getPeakGenes <- function(peak_anno) {
  if (is.null(peak_anno)) return(NULL)
  anno_df <- as.data.frame(peak_anno)
  gene_ids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
  return(gene_ids)
}

open_genes <- getPeakGenes(open_anno)
closed_genes <- getPeakGenes(closed_anno)
all_genes <- getPeakGenes(all_anno)

# GO 富集分析
runGOEnrichment <- function(genes, name) {
  if (is.null(genes) || length(genes) == 0) {
    cat(paste0("  跳过 ", name, " GO 分析 (无基因)\n"))
    return(NULL)
  }

  cat(paste0("  GO 分析: ", name, " (", length(genes), " genes)\n"))

  go_result <- enrichGO(
    gene = genes,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
    write.csv(as.data.frame(go_result), file.path(output_dir, paste0(name, "_GO_enrichment.csv")))

    # 柱状图
    p_bar <- barplot(go_result, showCategory = 20, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      ggtitle(paste0(name, " - GO Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_GO_barplot.png")), p_bar, width = 10, height = 12)

    # 气泡图
    p_dot <- dotplot(go_result, showCategory = 20, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      ggtitle(paste0(name, " - GO Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_GO_dotplot.png")), p_dot, width = 10, height = 12)
  }

  return(go_result)
}

go_open <- runGOEnrichment(open_genes, "open_peaks")
go_closed <- runGOEnrichment(closed_genes, "closed_peaks")
go_all <- runGOEnrichment(all_genes, "all_peaks")

# KEGG 富集分析
runKEGGEnrichment <- function(genes, name) {
  if (is.null(genes) || length(genes) == 0) {
    cat(paste0("  跳过 ", name, " KEGG 分析 (无基因)\n"))
    return(NULL)
  }

  cat(paste0("  KEGG 分析: ", name, "\n"))

  kegg_result <- enrichKEGG(
    gene = genes,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    kegg_result <- setReadable(kegg_result, OrgDb = org_db, keyType = "ENTREZID")

    write.csv(as.data.frame(kegg_result), file.path(output_dir, paste0(name, "_KEGG_enrichment.csv")))

    # 柱状图
    p_bar <- barplot(kegg_result, showCategory = 20) +
      ggtitle(paste0(name, " - KEGG Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_KEGG_barplot.png")), p_bar, width = 10, height = 8)

    # 气泡图
    p_dot <- dotplot(kegg_result, showCategory = 20) +
      ggtitle(paste0(name, " - KEGG Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_KEGG_dotplot.png")), p_dot, width = 10, height = 8)
  }

  return(kegg_result)
}

kegg_open <- runKEGGEnrichment(open_genes, "open_peaks")
kegg_closed <- runKEGGEnrichment(closed_genes, "closed_peaks")
kegg_all <- runKEGGEnrichment(all_genes, "all_peaks")

# ==================== Motif 分析准备 ====================

cat("\n准备 Motif 分析文件...\n")

# 导出 peak 序列（用于 MEME/HOMER motif 分析）
if (!is.null(open_gr)) {
  export.bed(open_gr, file.path(output_dir, "open_peaks_for_motif.bed"))
}
if (!is.null(closed_gr)) {
  export.bed(closed_gr, file.path(output_dir, "closed_peaks_for_motif.bed"))
}

# ==================== 汇总报告 ====================

cat("\n生成汇总报告...\n")

summary_df <- data.frame(
  Category = c("Open Peaks", "Closed Peaks", "Total Significant Peaks"),
  Count = c(nrow(open_peaks), nrow(closed_peaks), nrow(diff_peaks)),
  Genes_Annotated = c(length(open_genes), length(closed_genes), length(all_genes))
)

write.csv(summary_df, file.path(output_dir, "enrichment_summary.csv"), row.names = FALSE)

# 保存结果对象
save(open_anno, closed_anno, go_open, go_closed, kegg_open, kegg_closed,
     file = file.path(output_dir, "enrichment_results.RData"))

cat("\n==========================================\n")
cat("Peak 富集分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
