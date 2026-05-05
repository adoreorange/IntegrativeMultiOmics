# rawdata的数据质控和Harmony分析
rm(list=ls());gc()
setwd('/home/adore_org/B_scRNA-seq/analysis/')
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)

###### step1:导入数据 ######   
dir='RNA_seq_data/' 
fs=list.files('RNA_seq_data/',pattern = '*MWT')
fs


samples=list.files('RNA_seq_data/',pattern = '*MWT')
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1]
  pro2 = paste0(pro,'/filtered_feature_bc_matrix/')
  print(pro2)  
  tmp = Read10X(file.path(dir,pro2 )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  meta <- paste0('./barcode_info/', pro, '/Barcode_Infor.xls')
  metadata <- fread(meta, header = T , stringsAsFactors = F)
  rownames(metadata) <- metadata[,Barcode]
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  sce <- AddMetaData(object = sce, metadata = metadata)
  return(sce)
}) 
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples) 

sce.all[["RNA"]]@counts 

#看看合并前后的sce变化
sce.all <- readRDS('./2-harmony/sce.all_int.rds')
dim(sce.all[["RNA"]]@counts)


as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
# fivenum(sce.all$nFeature_RNA)
# table(sce.all$nFeature_RNA>800) 
# sce.all=sce.all[,sce.all$nFeature_RNA>800]
# sce.all

library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)

sp='mouse'
###### step2: QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 质量控制
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all)) # 原始数据的细胞数和基因数
print(dim(sce.all.filt)) # 质量控制后的细胞数和基因数
# 质量控制后的细胞减少了一点
setwd('../')
getwd()

# 查看质量控制后的细胞的ribo和mito比例
fivenum(sce.all.filt$percent_ribo)
fivenum(sce.all.filt$percent_mito)
table(sce.all.filt$nFeature_RNA> 500)

###### step3: harmony整合多个单细胞样品 ######
set.seed(10086)
table(sce.all$orig.ident)
table(sce.all.filt$orig.ident)
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all.filt)
  
  setwd('../')
  
}
sce <- sce.all.int
table(sce.all.filt$orig.ident)


# 查看高变基因
VariableFeatures(sce.all.int)

###### step4: 看标记基因库 ######
# 查看不同分辨率下的细胞分类
table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.2) 
table(sce.all.int$RNA_snn_res.0.8) 

getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.5')
setwd('check-by-0.5')
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.8')
setwd('check-by-0.8')
sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

last_markers_to_check
