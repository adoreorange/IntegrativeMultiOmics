#!/usr/bin/env Rscript
# 富集分析 (GO/KEGG)
# 用法: Rscript enrichment_analysis.R <deg_file> <output_dir>

# 加载必要的包
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)  # 人类基因注释，如果是其他物种请替换
  library(enrichplot)
  library(tidyverse)
  library(DOSE)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
deg_file <- ifelse(length(args) >= 1, args[1], "output/05_de_analysis/DEG_significant.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "output/06_enrichment")

# 可选：更换物种注释
# 模式生物注释包：
# 人类: org.Hs.eg.db
# 小鼠: org.Mm.eg.db
# 大鼠: org.Rn.eg.db
# 斑马鱼: org.Dr.eg.db
# 果蝇: org.Dm.eg.db
# 酵母: org.Sc.sgd.db

org_db <- org.Hs.eg.db

cat("==========================================\n")
cat("GO/KEGG 富集分析\n")
cat("==========================================\n")
cat("差异基因文件:", deg_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("==========================================\n")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取差异基因
cat("读取差异基因...\n")
deg <- read.csv(deg_file, row.names = 1)

# 准备基因列表
gene_list <- deg$log2FoldChange
names(gene_list) <- rownames(deg)
gene_list <- sort(gene_list, decreasing = TRUE)

# 上调基因和下调基因
up_genes <- names(gene_list[gene_list > 1])
down_genes <- names(gene_list[gene_list < -1])
all_deg_genes <- c(up_genes, down_genes)

cat("上调基因数:", length(up_genes), "\n")
cat("下调基因数:", length(down_genes), "\n")

# 基因 ID 转换（Symbol to Entrez）
cat("转换基因 ID...\n")
gene_df <- bitr(all_deg_genes, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org_db)
gene_entrez <- gene_df$ENTREZID

up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID

# ==================== GO 富集分析 ====================

cat("\n进行 GO 富集分析...\n")

# 所有 DEG 的 GO 分析
go_all <- enrichGO(
  gene = gene_entrez,
  OrgDb = org_db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 保存结果
if (nrow(as.data.frame(go_all)) > 0) {
  write.csv(as.data.frame(go_all), file = file.path(output_dir, "GO_enrichment_all.csv"))

  # GO 柱状图
  p_go_bar <- barplot(go_all, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    ggtitle("GO Enrichment Analysis")

  ggsave(file.path(output_dir, "GO_barplot.png"), p_go_bar, width = 10, height = 12)

  # GO 气泡图
  p_go_dot <- dotplot(go_all, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    ggtitle("GO Enrichment Analysis")

  ggsave(file.path(output_dir, "GO_dotplot.png"), p_go_dot, width = 10, height = 12)

  # GO 网络图
  p_go_cnet <- cnetplot(go_all, showCategory = 10, categorySize = "pvalue", foldChange = gene_list)
  ggsave(file.path(output_dir, "GO_cnetplot.png"), p_go_cnet, width = 12, height = 10)

  # GO 环形图
  p_go_circlize <- goplot(go_all, showCategory = 10)
  ggsave(file.path(output_dir, "GO_goplot.png"), p_go_circlize, width = 10, height = 10)
}

# 分别对上调和下调基因进行 GO 分析
if (length(up_entrez) > 0) {
  go_up <- enrichGO(gene = up_entrez, OrgDb = org_db, ont = "ALL", pvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go_up), file = file.path(output_dir, "GO_enrichment_up.csv"))
}

if (length(down_entrez) > 0) {
  go_down <- enrichGO(gene = down_entrez, OrgDb = org_db, ont = "ALL", pvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go_down), file = file.path(output_dir, "GO_enrichment_down.csv"))
}

# ==================== KEGG 富集分析 ====================

cat("\n进行 KEGG 富集分析...\n")

kegg_all <- enrichKEGG(
  gene = gene_entrez,
  organism = "hsa",  # 人类: hsa, 小鼠: mmu, 大鼠: rno
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (nrow(as.data.frame(kegg_all)) > 0) {
  # 转换基因 ID
  kegg_all <- setReadable(kegg_all, OrgDb = org_db, keyType = "ENTREZID")

  write.csv(as.data.frame(kegg_all), file = file.path(output_dir, "KEGG_enrichment_all.csv"))

  # KEGG 柱状图
  p_kegg_bar <- barplot(kegg_all, showCategory = 20) + ggtitle("KEGG Enrichment Analysis")
  ggsave(file.path(output_dir, "KEGG_barplot.png"), p_kegg_bar, width = 10, height = 8)

  # KEGG 气泡图
  p_kegg_dot <- dotplot(kegg_all, showCategory = 20) + ggtitle("KEGG Enrichment Analysis")
  ggsave(file.path(output_dir, "KEGG_dotplot.png"), p_kegg_dot, width = 10, height = 8)

  # KEGG 网络图
  p_kegg_cnet <- cnetplot(kegg_all, showCategory = 10, categorySize = "pvalue", foldChange = gene_list)
  ggsave(file.path(output_dir, "KEGG_cnetplot.png"), p_kegg_cnet, width = 12, height = 10)
}

# 分别对上调和下调基因进行 KEGG 分析
if (length(up_entrez) > 0) {
  kegg_up <- enrichKEGG(gene = up_entrez, organism = "hsa", pvalueCutoff = 0.05)
  if (nrow(as.data.frame(kegg_up)) > 0) {
    kegg_up <- setReadable(kegg_up, OrgDb = org_db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg_up), file = file.path(output_dir, "KEGG_enrichment_up.csv"))
  }
}

if (length(down_entrez) > 0) {
  kegg_down <- enrichKEGG(gene = down_entrez, organism = "hsa", pvalueCutoff = 0.05)
  if (nrow(as.data.frame(kegg_down)) > 0) {
    kegg_down <- setReadable(kegg_down, OrgDb = org_db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg_down), file = file.path(output_dir, "KEGG_enrichment_down.csv"))
  }
}

# ==================== GSEA 分析 ====================

cat("\n进行 GSEA 分析...\n")

# 转换基因列表
gene_list_entrez <- gene_list
names(gene_list_entrez) <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# GO GSEA
gsea_go <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org_db,
  ont = "ALL",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (nrow(as.data.frame(gsea_go)) > 0) {
  write.csv(as.data.frame(gsea_go), file = file.path(output_dir, "GSEA_GO_results.csv"))

  # GSEA 气泡图
  p_gsea <- dotplot(gsea_go, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    ggtitle("GSEA Analysis")
  ggsave(file.path(output_dir, "GSEA_dotplot.png"), p_gsea, width = 10, height = 12)

  # GSEA 山脊图
  p_ridge <- ridgeplot(gsea_go, showCategory = 15)
  ggsave(file.path(output_dir, "GSEA_ridgeplot.png"), p_ridge, width = 10, height = 8)
}

# 保存结果对象
save(go_all, kegg_all, gsea_go, file = file.path(output_dir, "enrichment_results.RData"))

cat("==========================================\n")
cat("富集分析完成!\n")
cat("结果保存在:", output_dir, "\n")
cat("==========================================\n")
