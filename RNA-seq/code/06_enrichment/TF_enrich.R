setwd("/home/adore_org/all_analysis/RNA/")
library(RcisTarget)
rm(list=ls());gc()

# Load gene sets to analyze. 
geneList_OF_YF <- read.csv('./rawdata/WT/OF_YF_WT.csv',header = T)
OF_YF_WT_up <- subset(geneList_OF_YF, padj < 0.05 & (log2FoldChange > 0.585))
OF_YF_WT_down <- subset(geneList_OF_YF, padj < 0.05 & (log2FoldChange < -0.585))
geneList_OM_YM <- read.csv('./rawdata/WT/OM_YM_WT.csv',header = T)
OM_YM_WT <- subset(geneList_OM_YM, padj < 0.05 & (log2FoldChange < -0.585 | log2FoldChange > 0.585))

setwd('./TF_enrich/')
geneLists <- list('OF_YF_up'=OF_YF_WT_up$symbol,
                  'OF_YF_down'=OF_YF_WT_down$symbol)


# Select motif database to use (i.e. organism and distance around TSS)
data("motifAnnotations_mgi_v9")
data(motifAnnotations_mgi)
motifRankings <- importRankings("/home/adore_org/scenic_ref/cisTarget_databases/mm10/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

# 合并进行
# Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_mgi_v9)

motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

# 分步进行
# motif 富集分析
# 1
motifs_AUC <- calcAUC(geneLists[1], motifRankings, nCores=15)

# 2
auc <- getAUC(motifs_AUC)["OF_YF_up",]
hist(auc, main="OF_YF_up", xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(auc)) + mean(auc)
abline(v=nes3, col="red")

# 3
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_mgi_v9)

class(motifEnrichmentTable)
dim(motifEnrichmentTable)

head(motifEnrichmentTable[,-"TF_lowConf", with=FALSE])

motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists)
dim(motifEnrichmentTable_wGenes)

# visules
geneSetName <- "OF_YF_up"
selectedMotifs <- c(sample(motifEnrichmentTable$motif, 2))
par(mfrow=c(2,2))
getSignificantGenes(geneLists[[geneSetName]], 
                    rankings=motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")


resultsSubset <- motifEnrichmentTable_wGenes[1:10,]
showLogo(resultsSubset)


anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

anotatedTfs$OF_YF_up


signifMotifNames <- motifEnrichmentTable$motif[1:4]

incidenceMatrix <- getSignificantGenes(geneLists$OF_YF_up, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
net <-visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                              nodesIdSelection = TRUE)

print(net)


library(visNetwork)
library(htmlwidgets)
library(webshot2)
saveWidget(net, "./TF_enrich/OF_YF_down_network.html")


# 将 HTML 转换为 PDF
webshot("./TF_enrich/network.html", "./TF_enrich/network.pdf", vwidth = 800, vheight = 600)
