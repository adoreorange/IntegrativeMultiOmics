# diffbind分析脚本：主要分析所有样本的差异Peak，包括PCAs、FRiP、均一化、对比分析、可视化
# all_diffbind.R
# input: bigWig
# output: diffbind result
setwd('/home/adore_org/diffbind/')
#BiocManager::install("lessR")
library(DiffBind, quietly = TRUE)
library(tidyverse, quietly = TRUE)
rm(list=ls());gc()

use_colors <- c(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

meta <- read.csv("./samples_meta.csv", header = TRUE)
##构建 dba 对象。
sampleDba <- dba(sampleSheet = meta)

#PCA
pdf('./plot/rawdata_pca.pdf',width = 6,height = 6)
plot(sampleDba)
dev.off()
#计算 count, 然后会给出 FRiP 数值。
sampleCount <- dba.count(sampleDba,bUseSummarizeOverlaps=T)
saveRDS(sampleCount,file=paste("./data_out/all/peak_data_collection.rds"))
sampleCount <- readRDS('./data_out/all/peak_data_collection.rds')

pdf('./plot/all/factor_pca.pdf',height = 9,width = 9)
dba.plotPCA(sampleCount,attributes = DBA_FACTOR,label = DBA_ID, vColors = use_colors,labelSize = 1)
dev.off()

# =======================================================
# show libsizes
info <- dba.show(sampleCount)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
write.csv(libsizes,'./data_out/all/libsizes.csv')

#进行均一化，然后分析差异的 Peaks.
sampleNor <- dba.normalize(sampleCount)
sampleNor
#计算 count, 然后会给出 FRiP 数值。
sampleCon <- dba.contrast(sampleNor,categories=DBA_TREATMENT,minMembers = 2)
#sampleCon <- dba.contrast(sampleNor,contrast = c("Factor","OF","YF"))
#sampleCon <- dba.contrast(sampleNor,contrast = c("Condition","Young","Old"))
#sampleCon <- dba.contrast(sampleNor, design = ~ Factor, contrast = c("Factor", "OM", "YM"), 
                          #reorderMeta = list(Factor="YM"))

con <- dba.show(DBA = sampleCon,bContrasts = T)
write.csv(con,'./data_out/all/contrasts.csv')
# test <-dba.contrast(tamoxifen,reorderMeta = list(Condition = "control"))
# dba.show(DBA = test,bContrasts = T)

sampleDiff <- dba.analyze(sampleCon,bBlacklist = FALSE, bGreylist = FALSE)
saveRDS(sampleDiff,'./data_out/all/all_sampleDiff.rds')
sampleDiff <- readRDS('./data_out/all/all_sampleDiff.rds')
dba.show(sampleDiff, bContrasts = TRUE)

plot(sampleDiff, contrast=5)
dba.plotMA(sampleDiff,contrast = 6)

saveRDS(sampleDiff, "./data_out/all/dbaAll.rds")
dbaAll <- sampleDiff
diffList <- list()
diffList$YF_OF <- dba.report(dbaAll, contrast = 1)
diffList$YF_YM <- dba.report(dbaAll, contrast = 2)
diffList$OF_OM <- dba.report(dbaAll, contrast = 5)
diffList$OM_YM <- dba.report(dbaAll, contrast = 6)

saveRDS(diffList, "./data_out/all/diffList.rds")

peakMeta <- as.data.table(dbaAll$peaks[[1]][, 1:3])

saveRDS(peakMeta, "./data_out/all/peakMeta.rds")

peakMtx <- map(dbaAll$peaks, ~ {.x$Score}) %>% reduce(cbind) %>% set_colnames(dbaAll$samples$SampleID) %>% as.data.table()

saveRDS(peakMtx, "./data_out/all/peakMtx.rds")



### YF_OF
diffPeaks_YF_OF <- dba.report(sampleDiff,th = 1,contrast = 1) %>%
  data.frame() %>%
  mutate(type = case_when(Fold >= 1 & FDR < 0.05 ~ "sigUp",
                          Fold <= -1 & FDR < 0.05 ~ "sigDown",
                          .default = "nonSig"))
# check
table(diffPeaks_YF_OF$type)
diffPeaks <- diffPeaks_YF_OF
# ==============================================================================
# annotate for peaks
# ==============================================================================
library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
peakAnno_all <- diffPeaks %>% GRanges()
peakAnno_all_plot <- annotatePeak(peakAnno_all, tssRegion = c(-3000, 3000),
                                  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)


peakAnno_all_df <- data.frame(peakAnno_all_plot)

# id transformation
ids <- mapIds(org.Mm.eg.db, keys = peakAnno_all_df$geneId, keytype = "ENTREZID", column = "SYMBOL")
peakAnno_all_df$GeneName <- ids
# check
table(peakAnno_all_df$type)

# export all diff
write.csv(peakAnno_all_df,file = paste('./data_out/OF_YF/diffbind_peaks_anno_diff_OF_YF.csv') ,row.names = F)

# plot 
peakAnno_df$annotation[grepl(x=peakAnno_df$annotation,pattern = 'Promoter')] <-'Promoter'

pdf(file = './plot/OF_YF/OF_YF_all_open_bar.pdf',height = 3,width = 12)
plotAnnoBar(peakAnno_all_plot)
dev.off()

pdf(file = './plot/OF_YF/OF_YF_all_open_pie.pdf',height = 6,width = 7)
plotAnnoPie(peakAnno_all_plot)
dev.off()

##### sig


sigUP <- diffPeaks %>% dplyr::filter(type == "sigUp") %>% GRanges()
sigDOWN <- diffPeaks %>% dplyr::filter(type == "sigDown") %>% GRanges()

# up ----
peakAnno_up_plot <- annotatePeak(sigUP, tssRegion = c(-3000, 3000),
                                 TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
peakAnno_up_df <- data.frame(peakAnno_up_plot)
ids <- mapIds(org.Mm.eg.db, keys = peakAnno_up_df$geneId, keytype = "ENTREZID", column = "SYMBOL")
peakAnno_up_df$GeneName <- ids
# output
write.csv(peakAnno_up_df, './data_out/OF_YF/OF_YF_sig_up.csv', row.names = FALSE)

# plot 启动子区域
peakAnno_up_df$annotation[grepl(x=peakAnno_up_df$annotation,pattern = 'Promoter')] <-'Promoter'

pdf(file = './plot/OF_YF/OF_YF_up_open_bar.pdf',height = 3,width = 12)
plotAnnoBar(peakAnno_up_plot)
dev.off()

pdf(file = './plot/OF_YF/OF_YF_up_open_pie.pdf',height = 6,width = 7)
plotAnnoPie(peakAnno_up_plot)
dev.off()
upgene <- data.frame(unique(peakAnno_up_df$GeneName))
write.csv(upgene,'./data_out/OF_YF/up_gene.csv')


# down ----
peakAnno_down_plot <- annotatePeak(sigDOWN, tssRegion = c(-3000, 3000),
                                   TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
peakAnno_down_df <- data.frame(peakAnno_down_plot)
ids <- mapIds(org.Mm.eg.db, keys = peakAnno_down_df$geneId, keytype = "ENTREZID", column = "SYMBOL")
peakAnno_down_df$GeneName <- ids
# output
write.csv(peakAnno_down_df, './data_out/OF_YF/OF_YF_sig_down.csv', row.names = FALSE)

# plot 启动子区域
peakAnno_down_df$annotation[grepl(x=peakAnno_down_df$annotation,pattern = 'Promoter')] <-'Promoter'

pdf(file = './plot/OF_YF/OF_YF_down_open_bar.pdf',height = 3,width = 12)
plotAnnoBar(peakAnno_down_plot)
dev.off()

pdf(file = './plot/OF_YF/OF_YF_down_open_pie.pdf',height = 6,width = 7)
plotAnnoPie(peakAnno_down_plot)
dev.off()
downgene <- data.frame(unique(peakAnno_down_df$GeneName))
write.csv(downgene,'./data_out/OF_YF/down_gene.csv')





