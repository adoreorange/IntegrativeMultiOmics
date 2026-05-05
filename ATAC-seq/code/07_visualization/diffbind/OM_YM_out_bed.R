# 输出OM_YM的差异Peak为bed文件

# ==============================================================================
# output peaks with types
# ==============================================================================
rm(list=ls());gc()
setwd('./out_bed/')
diffPeaks_OM_YM <- read.csv("../data_out/OM_YM/diffbind_peaks_anno_diff_OM_YM.csv")

diffPeaks <- diffPeaks_OM_YM
table(diffPeaks$type)
# check
types <- c("sigUp","sigDown","nonSig")

# x = 2
lapply(seq_along(types),function(x){
  tmp <- diffPeaks %>% dplyr::filter(type == types[x]) %>% dplyr::select(seqnames,start,end)
  
  write.table(tmp,file = paste0('./OM_YM/',types[x],".bed"),
              quote = F,col.names = F,row.names = F,sep = "\t")
})

'../bigwig/YF-WT.mRp.clN.bigWig'
# plot diff peaks region signal
nohup computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OM_YM/sigUp.bed ./out_bed/OM_YM/sigDown.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig_con/OM-WT_REP1.mLb.clN.bigWig ./bigwig_con/OM-WT_REP2.mLb.clN.bigWig ./bigwig_con/YM-WT_REP1.mLb.clN.bigWig ./bigwig_con/YM-WT_REP1.mLb.clN.bigWig \
--skipZeros -p 20 \
-o ./tss/OM_YM/OM_YM_matrix_diff.gz &

c('#E53831','#E43730','#DA3830','#E73B34')
plotHeatmap -m ./tss/OM_YM/OM_YM_matrix_diff.gz --zMin 0 --zMax 3 --yMin 0 --yMax 3.5 --colorList 'white,#0066CC' --heatmapHeight 16 --heatmapWidth 4 -o ./tss/OM_YM/OM_YM_diffPeaks.pdf


# plot diff peaks region signal
computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OF_YF_sigUp.bed ./out_bed/OF_YF_sigDown.bed ./out_bed/OF_YF_nonSig.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig/OF-WT.mRp.clN.bigWig ./bigwig/YF-WT.mRp.clN.bigWig \
--skipZeros -p 15 \
-o matrix_diff.gz