

# 1. Load packages --------------------------------------------------------
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/mycolors.R')
rm(list=ls());gc()
setwd("/home/adore_org/all_analysis/RNA/GSEA/")

library(fgsea)          
library(data.table)     
library(ggplot2)       
library(dplyr)          
library(msigdb)         
library(GSEABase) 
library(org.Mm.eg.db) 
library(clusterProfiler)
library(enrichplot)
library(GseaVis)



# 2. Load data ------------------------------------------------------------

# input data
OF_YF_WT <- read.csv('../Volcano/out_data/OF_WT_YF_WT_group.csv',row.names = 1)
# OF_YF_WT <- subset(OF_YF_WT, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))
OF_YF_WT <- OF_YF_WT[order(OF_YF_WT$log2FoldChange, decreasing = T), ]

table(OF_YF_WT$group)
head(OF_YF_WT)

# ID转化
gene_entrezid <- bitr(geneID = OF_YF_WT$symbol,
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", # 转成ENTREZID
                      OrgDb = "org.Mm.eg.db"
)
head(gene_entrezid)

gene_entrezid$logFC <- OF_YF_WT$log2FoldChange[match(gene_entrezid$SYMBOL, OF_YF_WT$symbol)]
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID 
head(genelist)

library(msigdbr)
gmt_file <- "./m5.all.v2024.1.Mm.entrez.gmt"  # 替换为你的GMT文件路径
gene_sets <- read.gmt(gmt_file)
head(gene_sets) # ENTREZID
gene_ids <- mapIds(org.Mm.eg.db, keys = gene_sets$gene, column = "SYMBOL", keytype = "ENTREZID")
gene_sets$gene <- gene_ids

#database <- msigdbr(species= "Mus musculus") %>% dplyr::select(gs_name, entrez_gene)
#head(database)

# GSEA GO
res <- GSEA(genelist, 
            TERM2GENE = gene_sets,
            minGSSize = 5, maxGSSize = 500,
            pvalueCutoff = 1, pAdjustMethod = "BH")

res<- setReadable(res, OrgDb= org.Mm.eg.db, keyType= "ENTREZID")
result <- as.data.frame(gsea_res)
result <-subset(result,pvalue<0.05)
result <- result[order(result$NES, decreasing = T), ]
write.csv(result,'./OF_WT_YF_WT_GSEA_GO.csv')

# GSEA KEGG
KEGG_ges<- gseKEGG(geneList= genelist, organism= "mmu", 
                   minGSSize= 5, maxGSSize= 500, 
                   pvalueCutoff= 1, pAdjustMethod= "BH", 
                   verbose= FALSE)

KEGG_ges<- setReadable(KEGG_ges, OrgDb= org.Mm.eg.db, keyType= "ENTREZID")
result_KEGG <- as.data.frame(KEGG_ges)
result_KEGG <-subset(result_KEGG,pvalue<0.05)
result_KEGG <- result_KEGG[order(result_KEGG$NES, decreasing = T), ]
write.csv(result_KEGG,'./OF_WT_YF_WT_GSEA_KEGG.csv')

# plot GSEA GO
source('./GSEAplot.R')
res = gsea_res
p <- plot_gseplot(data = NULL, data_ud = result, x = 11);p
ggsave("./plot/GOBP_RESPONSE_TO_INTERFERON_BETA.pdf", p, height = 5,width = 8)

# plot GSEA KEGG
res = KEGG_ges
p <- plot_gseplot(data = NULL, data_ud = result_KEGG, x = 4);p
ggsave("./plot/PI3K-Akt signaling pathway.pdf", p, height = 5,width = 8)



# mutliplot 
gseaplot2(gsea_res, geneSetID = rownames(result)[1:4],
          color = mycolors[1:4], base_size = 12, rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3, pvalue_table =  F, ES_geom = 'line')


# 山脊图
library(enrichplot)
library(ggplot2)

ridgeplot(gsea_res,
          showCategory = 20,
          fill = "pvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30, orderBy = "NES", decreasing = FALSE) + 
  theme(axis.text.y = element_text(size=8))

ids <- gsea_res@result$ID[1:5]

gseadist(gsea_res,
         IDs = ids,
         type = "density" # boxplot
) + theme(legend.direction = "vertical")


gseaplot(gsea_res, geneSetID = 2, by = "runningScore", 
         title = gsea_res$Description[2])

gseaplot(gsea_res, geneSetID = 1, by = "preranked", 
         title = gsea_res$Description[1]) + 
  theme(plot.title = element_text(size = 10, color = "blue"))

# 如果两个子图都画的话返回的是一个ggplist对象，此时如果要修改图形细节，可以使用取子集的方法提取其中的子图形，此时的子图形是ggplot对象，又可以使用ggplot2语法修改了
p <- gseaplot(gsea_res, geneSetID = 1, title = gsea_res$Description[1])
p

#取子集进行修改
p[[1]] <- p[[1]]+theme(plot.title = element_text(size = 6))
p



# 默认subplots = 1:3，把3个图放一起
gseaplot2(gsea_res,geneSetID = 1,title = "title",
          subplots = 1:3,
          base_size = 10)

gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
gseaplot2(gsea_res, geneSetID = 1, subplots = 1:2)

#把entrezid变为symbol
gsea_res_symbol <- setReadable(gsea_res, "org.Mm.eg.db", "ENTREZID")

p <- gseaplot2(gsea_res_symbol,geneSetID = 6, title = gsea_res_symbol$Description[6])

p[[1]] <- p[[1]]+ theme(title = element_text(color = "red"))
p
