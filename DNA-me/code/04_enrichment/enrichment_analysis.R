#!/usr/bin/env Rscript
# DMR 功能富集分析
# 用法: Rscript enrichment_analysis.R <dmr_file> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(annotatr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ChIPseeker)
  library(tidyverse)
  library(ggplot2)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
dmr_file <- ifelse(length(args) >= 1, args[1], "output/03_dmr_analysis/DMR_significant.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "output/04_enrichment")

# 物种设置
org_db <- org.Hs.eg.db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- "hg38"

cat("==========================================\n")
cat("DMR 功能富集分析\n")
cat("==========================================\n")
cat("DMR 文件:", dmr_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==================== 读取 DMR ====================

cat("\n读取 DMR 数据...\n")
dmr_data <- read.csv(dmr_file, stringsAsFactors = FALSE)

cat("DMR 数量:", nrow(dmr_data), "\n")

# 创建 GRanges
dmr_gr <- GRanges(
  seqnames = dmr_data$chr,
  ranges = IRanges(start = dmr_data$start, end = dmr_data$end),
  direction = dmr_data$direction,
  mean_diff = dmr_data$mean_diff
)

# 分别处理高甲基化和低甲基化区域
hyper_gr <- dmr_gr[dmr_gr$direction == "Hyper"]
hypo_gr <- dmr_gr[dmr_gr$direction == "Hypo"]

cat("高甲基化 DMR:", length(hyper_gr), "\n")
cat("低甲基化 DMR:", length(hypo_gr), "\n")

# ==================== 基因组特征注释 ====================

cat("\n进行基因组特征注释...\n")

# 使用 annotatr 进行注释
# 构建注释
annots <- c(
  paste0(genome, "_genes_cpgs"),
  paste0(genome, "_genes_exons"),
  paste0(genome, "_genes_firstexons"),
  paste0(genome, "_genes_intronsexons"),
  paste0(genome, "_genes_introns"),
  paste0(genome, "_genes_promoters"),
  paste0(genome, "_cpgs")
)

# 注释函数
annotateDMR <- function(gr, name) {
  if (length(gr) == 0) {
    cat(paste0("  跳过 ", name, " (无 DMR)\n"))
    return(NULL)
  }

  cat(paste0("  注释: ", name, "\n"))

  # 使用 ChIPseeker 注释
  peak_anno <- annotatePeak(
    gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )

  # 保存注释结果
  anno_df <- as.data.frame(peak_anno)
  write.csv(anno_df, file.path(output_dir, paste0(name, "_annotation.csv")), row.names = FALSE)

  # 基因组特征分布图
  p_anno <- plotAnnoBar(peak_anno) +
    ggtitle(paste0(name, " - Genomic Annotation")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, paste0(name, "_genomic_annotation.png")), p_anno, width = 10, height = 6)

  # TSS 距离分布
  p_dist <- plotDistToTSS(peak_anno) +
    ggtitle(paste0(name, " - Distance to TSS")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, paste0(name, "_distance_to_tss.png")), p_dist, width = 10, height = 6)

  return(peak_anno)
}

# 注释所有 DMR
hyper_anno <- annotateDMR(hyper_gr, "hyper_DMR")
hypo_anno <- annotateDMR(hypo_gr, "hypo_DMR")
all_anno <- annotateDMR(dmr_gr, "all_DMR")

# ==================== CpG 岛分析 ====================

cat("\n分析 CpG 岛分布...\n")

# 获取 CpG 岛注释
cpg_islands <- getGenomicAnnotation(gr = dmr_gr, TxDb = txdb, annoDb = "org.Hs.eg.db")

# CpG 上下文统计
cpg_context <- data.frame(
  Region = c("CpG Island", "Shore", "Shelf", "Open Sea"),
  Description = c(
    "CpG islands (CGI)",
    "0-2kb from CGI",
    "2-4kb from CGI",
    ">4kb from CGI"
  )
)

write.csv(cpg_context, file.path(output_dir, "cpg_context_info.csv"), row.names = FALSE)

# ==================== 基因富集分析 ====================

cat("\n进行基因富集分析...\n")

# 提取基因 ID
getGeneIds <- function(peak_anno) {
  if (is.null(peak_anno)) return(NULL)
  anno_df <- as.data.frame(peak_anno)
  gene_ids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
  return(gene_ids)
}

hyper_genes <- getGeneIds(hyper_anno)
hypo_genes <- getGeneIds(hypo_anno)
all_genes <- getGeneIds(all_anno)

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

go_hyper <- runGOEnrichment(hyper_genes, "hyper_DMR")
go_hypo <- runGOEnrichment(hypo_genes, "hypo_DMR")
go_all <- runGOEnrichment(all_genes, "all_DMR")

# KEGG 富集分析
runKEGGEnrichment <- function(genes, name) {
  if (is.null(genes) || length(genes) == 0) {
    cat(paste0("  跳过 ", name, " KEGG 分析\n"))
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

    # 可视化
    p_bar <- barplot(kegg_result, showCategory = 20) +
      ggtitle(paste0(name, " - KEGG Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_KEGG_barplot.png")), p_bar, width = 10, height = 8)

    p_dot <- dotplot(kegg_result, showCategory = 20) +
      ggtitle(paste0(name, " - KEGG Enrichment"))

    ggsave(file.path(output_dir, paste0(name, "_KEGG_dotplot.png")), p_dot, width = 10, height = 8)
  }

  return(kegg_result)
}

kegg_hyper <- runKEGGEnrichment(hyper_genes, "hyper_DMR")
kegg_hypo <- runKEGGEnrichment(hypo_genes, "hypo_DMR")
kegg_all <- runKEGGEnrichment(all_genes, "all_DMR")

# ==================== 通路分析 ====================

cat("\n进行通路分析...\n")

# Reactome 富集（如果安装了 ReactomePA）
if (requireNamespace("ReactomePA", quietly = TRUE)) {
  library(ReactomePA)

  runReactomeEnrichment <- function(genes, name) {
    if (is.null(genes) || length(genes) == 0) return(NULL)

    cat(paste0("  Reactome 分析: ", name, "\n"))

    reactome_result <- enrichPathway(
      gene = genes,
      organism = "human",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )

    if (!is.null(reactome_result) && nrow(as.data.frame(reactome_result)) > 0) {
      write.csv(as.data.frame(reactome_result), file.path(output_dir, paste0(name, "_Reactome_enrichment.csv")))

      p_dot <- dotplot(reactome_result, showCategory = 20) +
        ggtitle(paste0(name, " - Reactome Pathway"))

      ggsave(file.path(output_dir, paste0(name, "_Reactome_dotplot.png")), p_dot, width = 10, height = 8)
    }

    return(reactome_result)
  }

  reactome_all <- runReactomeEnrichment(all_genes, "all_DMR")
}

# ==================== 汇总报告 ====================

cat("\n生成汇总报告...\n")

summary_df <- data.frame(
  Category = c("Hyper-methylated DMR", "Hypo-methylated DMR", "Total DMR"),
  Count = c(length(hyper_gr), length(hypo_gr), length(dmr_gr)),
  Annotated_Genes = c(length(hyper_genes), length(hypo_genes), length(all_genes))
)

write.csv(summary_df, file.path(output_dir, "enrichment_summary.csv"), row.names = FALSE)

# 保存结果对象
save(hyper_anno, hypo_anno, go_hyper, go_hypo, kegg_hyper, kegg_hypo,
     file = file.path(output_dir, "enrichment_results.RData"))

cat("\n==========================================\n")
cat "富集分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
