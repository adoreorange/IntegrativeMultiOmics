# 输出所有样本的差异Peak为bed文件

# ==============================================================================
# output peaks with types
# ==============================================================================
setwd('./out_bed/')
diffPeaks <- read.csv("../data_out/diffbind_peaks_diff_OF_YF.csv")

types <- c("sigUp","sigDown","nonSig")

# x = 2
lapply(seq_along(types),function(x){
  tmp <- diffPeaks %>% dplyr::filter(type == types[x]) %>% dplyr::select(seqnames,start,end)
  
  write.table(tmp,file = paste0('./out_bed/OF_YF/',types[x],".bed"),
              quote = F,col.names = F,row.names = F,sep = "\t")
})

'../bigwig/YF-WT.mRp.clN.bigWig'
# plot diff peaks region signal
computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OF_YF_sigUp.bed ./out_bed/OF_YF_sigDown.bed ./out_bed/OF_YF_nonSig.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig/OF-WT.mRp.clN.bigWig ./bigwig/YF-WT.mRp.clN.bigWig \
--skipZeros -p 15 \
-o matrix_diff.gz

plotHeatmap -m matrix_diff.gz --colorList 'white,#0066CC' --heatmapHeight 12 --heatmapWidth 6 -o diffPeaks.pdf


# plot diff peaks region signal
computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OF_YF_sigUp.bed ./out_bed/OF_YF_sigDown.bed ./out_bed/OF_YF_nonSig.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig/OF-WT.mRp.clN.bigWig ./bigwig/YF-WT.mRp.clN.bigWig \
--skipZeros -p 15 \
-o matrix_diff.gz