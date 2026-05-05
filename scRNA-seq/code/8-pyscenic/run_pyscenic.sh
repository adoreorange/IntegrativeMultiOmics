# 运行pyscenic分析
CURRENT_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo "The current time is: $CURRENT_TIME"
echo '# 开始分析'
# 不同物种的数据库不一样，这里是小鼠是 mouse 
dir=/home/huyifeng/pyscenic #改成自己的目录
tfs=$dir/TFs/mm9/allTFs_mm.txt
feather=$dir/TFs/mm9/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
feather2=$dir/TFs/mm9/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/TFs/mm9/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl 
# 一定要保证上面的数据库文件完整无误哦 
input_loom=$dir/data/Breg.loom
ls $tfs  $feather $feather2  $tbl  

# pyscenic 的3个步骤之 grn
echo '步骤一 grn'
pyscenic grn \
	--num_workers 20 \
	--output adj.sample.tsv \
	--method grnboost2 \
	$input_loom  \
	$tfs 

#pyscenic 的3个步骤之 cistarget
echo '步骤二 cistarget'
pyscenic ctx \
	adj.sample.tsv $feather2 \
	--annotations_fname $tbl \
	--expression_mtx_fname $input_loom  \
	--mode "dask_multiprocessing" \
	--output Breg.csv \
	--num_workers 20  \
	--mask_dropouts

#pyscenic 的3个步骤之 AUCell
echo '步骤三 AUcell'
pyscenic aucell \
	$input_loom \
	Breg.csv \
	--output out_SCENIC.loom \
	--num_workers 20

CURRENT_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo "处理完毕"
echo "The current time is: $CURRENT_TIME"
