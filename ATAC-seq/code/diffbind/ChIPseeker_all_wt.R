setwd('/home/adore_org/diffbind/')
# annotate for peaks
# ==============================================================================
library(ChIPseeker)
library(clusterProfiler)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

## 导入文件
OF1 <- readPeakFile("./beds/OF-WT_REP1.mLb.clN_summits.bed")
OF2 <- readPeakFile("./beds/OF-WT_REP2.mLb.clN_summits.bed")
YF1 <- readPeakFile("./beds/YF-WT_REP1.mLb.clN_summits.bed")
YF2 <- readPeakFile("./beds/YF-WT_REP2.mLb.clN_summits.bed")
OM1 <- readPeakFile("./beds/OM-WT_REP1.mLb.clN_summits.bed")
OM2 <- readPeakFile("./beds/OM-WT_REP2.mLb.clN_summits.bed")
YM1 <- readPeakFile("./beds/YM-WT_REP1.mLb.clN_summits.bed")
YM2 <- readPeakFile("./beds/YM-WT_REP2.mLb.clN_summits.bed")
peaks <- list(OF1=OF1,OF2=OF2,YF1=YF1,YF2=YF2,
              OM1=OM1,OM2=OM2,YM1=YM1,YM2=YM2)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,
                       addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")

pdf('./plot/OF_YF/Feature Distribution_OF_YF_WT.pdf', width = 12, height = 8, useDingbats = T)
plotAnnoBar(peakAnnoList[1:4])
dev.off()

pdf('./plot/OF_YF/Feature_Distribution_TSS_OF_YF_WT.pdf', width = 12, height = 8, useDingbats = T)
plotDistToTSS(peakAnnoList[1:4],title="Feature Distribution relative to TSS")
dev.off()

pdf('./plot/OM_YM/Feature Distribution_OM_YM_WT.pdf', width = 12, height = 8, useDingbats = T)
plotAnnoBar(peakAnnoList[5:8])
dev.off()

pdf('./plot/OM_YM/Feature_Distribution_TSS_OM_YM_WT.pdf', width = 12, height = 8, useDingbats = T)
plotDistToTSS(peakAnnoList[5:8],title="Feature Distribution relative to TSS")
dev.off()

# Output peakAnnolist file
save(peakAnnoList,file="./data_out/All_WT_peakAnnolist.rda")
saveRDS(peakAnnoList,file="./data_out/All_WT_peakAnnolist.rds")



write.table(as.data.frame(peakAnnoList$OF1),file="./data_out/All_WT.PeakAnno",sep='	',quote = F)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Nanog.csv")
