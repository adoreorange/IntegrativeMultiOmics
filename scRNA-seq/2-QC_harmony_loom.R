rm(list=ls());gc()
setwd('/home/adore_org/B_scRNA-seq/integrated//')
options(stringsAsFactors = F) 
source('./scRNA_scripts/lib.R')
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)

###### step1:导入数据 ######   
# before filter
RNA <- readRDS('./pydata/luo_loom.rds')
dim(RNA[["RNA"]]@counts )


dir='GSE155006_RAW/'
fs=list.files('./GSE155006_RAW/',pattern = '*')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]
samples
##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
if(F){
lapply(unique(samples),function(x){
  #x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE155006_RAW/", paste(str_split(y[1],'_',simplify = T)[,4], collapse = "_"))
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE155006_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE155006_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE155006_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
}

dir='GSE155006_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
}) 
samples <-c('Aged-17','Young-3')
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  ) 



sp='mouse'
###### step2: QC质控 ######
dir.create("./1-loom-QC")
setwd("./1-loom-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
RNA.all.filt = basic_qc(RNA)
print(dim(RNA));print(dim(RNA.all.filt))
setwd('../')
getwd()

###### step2: QC质控 ######
dir.create("./1-QC-GSE")
setwd("./1-QC-GSE")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()

fivenum(sce.all.filt$percent_ribo)
table(sce.all.filt$nFeature_RNA> 5)


###### step3: harmony整合多个单细胞样品 ######
set.seed(10086)
table(RNA.all.filt$orig.ident)
if(T){
  dir.create("2-harmony_loom")
  getwd()
  setwd("2-harmony_loom")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  
  RNA.all.int = run_harmony(RNA)
  setwd('../')
  getwd()
}

# 查看高变基因
VariableFeatures(RNA.all.int)
sce.all.filt$Sample <- sce.all.filt$orig.ident
sce.all.filt$order <- substr(x = colnames(sce.all.filt),start = 1,stop = 7)

con.all=merge(x=sce.all.filt,
              y=RNA.all.filt) 
