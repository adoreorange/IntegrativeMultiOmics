# footprinting分析脚本
# input: sorted.bam, narrowPeak
# output: footprinting result

date
cd /home/huyifeng/ATAC/con_rep/
pwd
for id in OF-WT OM-WT YF-WT YM-WT
do
	echo $id
	rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=./ --output-prefix=$id \
	       	${id}.mRp.clN.sorted.bam ${id}.mRp.clN_peaks.narrowPeak
done
echo -e " \n \n \n  333# ALL  Work Done !!!\n \n \n "
date
