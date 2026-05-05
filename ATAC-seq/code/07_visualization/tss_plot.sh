# tss plot脚本
# input: bigWig
# output: tss plot

nohup computeMatrix reference-point --referencePoint TSS -R genes_Gencode_mm10.bed -b 3000 -a 3000 --missingDataAsZero -S ./bigwig/OF-WT_REP1.mLb.clN.bigWig ./bigwig/OF-WT_REP2.mLb.clN.bigWig ./bigwig/YF-WT_REP1.mLb.clN.bigWig ./bigwig/YF-WT_REP2.mLb.clN.bigWig ./bigwig/OM-WT_REP1.mLb.clN.bigWig ./bigwig/OM-WT_REP2.mLb.clN.bigWig ./bigwig/YM-WT_REP1.mLb.clN.bigWig ./bigwig/YM-WT_REP2.mLb.clN.bigWig  ./bigwig/OF-KO_REP1.mLb.clN.bigWig ./bigwig/OF-KO_REP2.mLb.clN.bigWig ./bigwig/YF-KO_REP1.mLb.clN.bigWig ./bigwig/YF-KO_REP2.mLb.clN.bigWig ./bigwig/OM-KO_REP1.mLb.clN.bigWig ./bigwig/OM-KO_REP2.mLb.clN.bigWig ./bigwig/YM-KO_REP1.mLb.clN.bigWig ./bigwig/YM-KO_REP2.mLb.clN.bigWig --skipZeros -p 20 -o ./com_bw/all_3kbp.gz &