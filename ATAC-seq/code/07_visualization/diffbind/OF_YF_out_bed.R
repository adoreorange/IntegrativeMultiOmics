# 输出OF_YF的差异Peak为bed文件

# ==============================================================================
# output peaks with types
# ==============================================================================
setwd('./out_bed/')
diffPeaks_OF_YF <- read.csv("../data_out/OF_YF/diffbind_peaks_anno_diff_OF_YF.csv")

diffPeaks <- diffPeaks_OF_YF
table(diffPeaks$type)

# check
types <- c("sigUp","sigDown","nonSig")

# x = 2
lapply(seq_along(types),function(x){
  tmp <- diffPeaks %>% dplyr::filter(type == types[x]) %>% dplyr::select(seqnames,start,end)
  
  write.table(tmp,file = paste0('./OF_YF/',types[x],".bed"),
              quote = F,col.names = F,row.names = F,sep = "\t")
})

'../bigwig/YF-WT.mRp.clN.bigWig'
# plot diff peaks region signal
nohup computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OF_YF/sigUp.bed ./out_bed/OF_YF/sigDown.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig_con/OF-WT_REP1.mLb.clN.bigWig ./bigwig_con/OF-WT_REP2.mLb.clN.bigWig ./bigwig_con/YF-WT_REP1.mLb.clN.bigWig ./bigwig_con/YF-WT_REP1.mLb.clN.bigWig \
--skipZeros -p 20 \
-o ./tss/OF_YF/OF_YF_matrix_diff.gz &

c('#E53831','#E43730','#DA3830','#E73B34')
plotHeatmap -m ./tss/OF_YF/OF_YF_matrix_diff.gz --zMin 0 --zMax 3 --yMin 0 --yMax 3.5 --colorList 'white,#E43730' --heatmapHeight 16 --heatmapWidth 4 -o ./tss/OF_YF/OF_YF_diffPeaks.pdf

findMotifsGenome.pl  /home/diffbind/out_bed/OF_YF/sigUp.bed  mm10  ./OF_YF_OUT  -len 8,10,12

# plot diff peaks region signal
computeMatrix reference-point --referencePoint TSS \
-R ./out_bed/OF_YF_sigUp.bed ./out_bed/OF_YF_sigDown.bed ./out_bed/OF_YF_nonSig.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-S ./bigwig/OF-WT.mRp.clN.bigWig ./bigwig/YF-WT.mRp.clN.bigWig \
--skipZeros -p 15 \
-o matrix_diff.gz