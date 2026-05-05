# diffbind分析脚本：主要分析所有样本的差异Peak，包括PCAs、FRiP、均一化、对比分析、可视化
# OF_YF.R
# input: bigWig
# output: diffbind result
# 分析OF_YF的差异Peak
setwd('/home/adore_org/diffbind/')
#BiocManager::install("profileplyr")
library(DiffBind, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(profileplyr)
rm(list=ls());gc()

use_colors <- c(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

meta <- read.csv("./samples_meta.csv", header = TRUE)
##构建 dba 对象。
sampleDba <- dba(sampleSheet = meta)

#PCA
pdf('./plot/OF_YF/rawdata_pca.pdf',width = 6,height = 6)
plot(sampleDba)
dev.off()
#计算 count, 然后会给出 FRiP 数值。
sampleCount <- dba.count(sampleDba,bUseSummarizeOverlaps = T)
saveRDS(sampleCount,file="./data_out/OF_YF/OF_YF_sampleCount.rds")
sampleCount <- readRDS('./data_out/OF_YF/OF_YF_sampleCount.rds')

pdf('./plot/factor_pca.pdf',height = 9,width = 9)
dba.plotPCA(sampleCount,attributes = DBA_FACTOR,label = DBA_ID, vColors = use_colors,labelSize = 1)
dev.off()

# =======================================================
# show libsizes
info <- dba.show(sampleCount)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
write.csv(libsizes,'./data_out/OF_YF/OF_YF_libsizes.csv')

#进行均一化，然后分析差异的 Peaks.
sampleNor <- dba.normalize(sampleCount)
#计算 count, 然后会给出 FRiP 数值。
#sampleCon <- dba.contrast(sampleNor,contrast = c("Factor","OF","YF"))
#sampleCon <- dba.contrast(sampleNor,contrast = c("Condition","Young","Old"))
sampleCon <- dba.contrast(sampleNor, design = ~ Factor, contrast = c("Factor", "OF", "YF"), 
                          reorderMeta = list(Factor="YF"))
dba.show(DBA = sampleCon,bContrasts = T)

# test <-dba.contrast(tamoxifen,reorderMeta = list(Condition = "control"))
# dba.show(DBA = test,bContrasts = T)

sampleDiff <- dba.analyze(sampleCon, bBlacklist = FALSE, bGreylist = FALSE)
dba.show(sampleDiff, bContrasts = TRUE)
saveRDS(sampleDiff,'./data_out/OF_YF/OF_YF_sampleDiff.rds')

sampleDiff <- readRDS('./data_out/OF_YF/OF_YF_sampleDiff.rds')


diffPeaks <- dba.report(sampleDiff,th = 1) %>%
  data.frame() %>%
  mutate(type = case_when(Fold >= 2 & FDR < 0.05 ~ "sigUp",
                          Fold <= -1 & FDR < 0.05 ~ "sigDown",
                          .default = "nonSig"))

# check
table(diffPeaks$type)
# ==============================================================================
# annotate for peaks
# ==============================================================================
library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
peakAnno_all <- diffPeaks %>% GRanges()
peakAnno_all_plot <- annotatePeak(peakAnno_all, tssRegion = c(-3000, 3000),assignGenomicAnnotation = T,
                                  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

plotAnnoBar(peakAnno_all_plot)
plotDistToTSS(peakAnno_all_plot)
vennpie(peakAnno_all_plot)
plotAnnoPie(peakAnno_all_plot)

peakAnno_all_df <- data.frame(peakAnno_all_plot)

# id transformation
SYMBOL <- mapIds(org.Mm.eg.db, keys = peakAnno_all_df$geneId, keytype = "ENTREZID", column = "SYMBOL")
GENENAME <- mapIds(org.Mm.eg.db, keys = peakAnno_all_df$geneId, keytype = "ENTREZID", column = "GENENAME")
peakAnno_all_df$SYMBOL <- SYMBOL
peakAnno_all_df$GENENAME <- GENENAME
# check
table(peakAnno_all_df$type)

# export all diff
write.csv(peakAnno_all_df,file = './data_out/OF_YF/diffbind_peaks_anno_diff_OF_YF.csv' ,row.names = F)

# plot

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

# plot
# peakAnno_up_df$annotation[grepl(x=peakAnno_up_df$annotation,pattern = 'Promoter')] <-'Promoter'

pdf(file = './plot/OF_YF/OF_YF_up_open_bar.pdf',height = 3,width = 12)
plotAnnoBar(peakAnno_up_plot)
dev.off()

pdf(file = './plot/OF_YF/OF_YF_up_open_pie.pdf',height = 6,width = 7)
plotAnnoPie(peakAnno_up_plot)
dev.off()
upgene <- data.frame(unique(peakAnno_up_df$GeneName))
write.csv(upgene,'./data_out/OF_YF/OF_YF_up_gene.csv')
plotAnnoBar(peakAnno_up_plot)

# down ----
peakAnno_down_plot <- annotatePeak(sigDOWN, tssRegion = c(-3000, 3000),
                                   TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
peakAnno_down_df <- data.frame(peakAnno_down_plot)
ids <- mapIds(org.Mm.eg.db, keys = peakAnno_down_df$geneId, keytype = "ENTREZID", column = "SYMBOL")
peakAnno_down_df$GeneName <- ids
# output
write.csv(peakAnno_down_df, './data_out/OF_YF/OF_YF_sig_down.csv', row.names = FALSE)

# plot
# peakAnno_down_df$annotation[grepl(x=peakAnno_down_df$annotation,pattern = 'Promoter')] <-'Promoter'

pdf(file = './plot/OF_YF/OF_YF_down_open_bar.pdf',height = 3,width = 12)
plotAnnoBar(peakAnno_down_plot)
dev.off()

pdf(file = './plot/OF_YF/OF_YF_down_open_pie.pdf',height = 6,width = 7)
plotAnnoPie(peakAnno_down_plot)
dev.off()
downgene <- data.frame(unique(peakAnno_down_df$GeneName))
write.csv(downgene,'./data_out/OF_YF/OF_YF_down_gene.csv')





