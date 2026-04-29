echo $date
input_loom=/home/huyifeng/pyscenic/Breg.loom
reg_csv=/home/huyifeng/pyscenic/reg.csv
#pyscenic 的3个步骤之 AUCell
pyscenic aucell \
	$input_loom \
	$reg_csv \
	--output out_SCENIC.loom \
	--num_workers 20 
echo $date
