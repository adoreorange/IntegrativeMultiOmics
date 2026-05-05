# 划分B细胞群，区分B1和B2细胞

rm(list=ls());gc()
setwd('/home/adore_org/B_scRNA-seq/analysis/')
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')

# creat dir
dir.create("4-B")
setwd("4-B/") 
getwd()

# creat data
set.seed(12345)
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
table(scRNA$celltype)
sce1 = scRNA[, scRNA$celltype %in% c( 'B' )]
cellinfo <- subset(sce1@meta.data, select = c("orig.ident", "percent_mito", "Sample", "Barcode",'percent_ribo','percent.redcell'))
sce <- CreateSeuratObject(sce1@assays$RNA@counts, meta.data = cellinfo)

as.data.frame(sce@assays$RNA@counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident) 

# harmony data

set.seed(10086)
table(sce$orig.ident)
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.int = run_harmony(sce)
  setwd('../')
  
}
sce <- sce.int
sce=readRDS( "./2-harmony/sce.all_int.rds")
sel.clust = "RNA_snn_res.0.5"
sce <- SetIdent(sce, value = sel.clust)

# plot res_tree
p2_tree=clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 9,width = 9)

# plot dimplot
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'umap',group.by = 'RNA_snn_res.0.5');p
ggsave(plot=p, filename="B_cell_umap_res0.5.pdf",height = 6,width = 7)
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'tsne',group.by = 'RNA_snn_res.0.5');p
ggsave(plot=p, filename="B_cell_tsne_res0.5.pdf",height = 6,width = 7)

### 根据基因marker划分B1和B2细胞
genes_to_check = c('Zbtb32','Bhlhe41','Fcrl5','Zcwpw1','Cd5', # B1 cells
                   'Ighd','Fcer2a','Ms4a4c','Vpreb3','Cr2','Cd1d1') #B2 ells
library(stringr)  
genes_to_check=str_to_title(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),cols = c("#FFFF00","#FF0000"),
            assay='RNA'  )  + coord_flip();p
ggsave('Breg_gene_markers.pdf',height = 10,width = 7)

##左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)

axes <- result$axes
label <- result$label

umap =DimPlot(sce, reduction = "umap",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.3",label = T,label.box = T,split.by = 'Sample') +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank());umap
ggsave('RNA_snn_res.0.3_umap_sample.pdf',width = 18,height = 7)

# canlcate marker genes
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/myfindmarkers.R')
dir.create('sample')
setwd('sample')
Idents(sce) <- sce$RNA_snn_res.0.5
myfindmarkers(RNA,filename = 'sample',min.pct = 0.1,species = 'mouse',colors = my36colors,)
setwd('../')
getwd()

# stacked_violin_plot
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/stacked_violin_plot.R')
stacked_violin_plot(gene = genes_to_check,seurat_object = sce,text.size = 10,flip = F,col = mycolors,
                    filename = "B1_B2_vlnplot",width = 12,height = 10,limits.max = 9,Mean=F)



# rationplot
source('/home/adore_org/Breg/scRNA_scripts/myrationplot.R')
myrationplot(seurat_object = sce,filename = 'B',col = mycolors)

# dotplot
p = DotPlot(sce, features=genes_to_check, assay = 'RNA',cols = mycolors) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,size = 10),
        axis.text.y = element_text(size = 12))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expressed") )+ #legend
  scale_color_gradientn(colours = c("#FFFF00","#FF0000"));p
#ggsave('check_B_markers.pdf',height = 11,width = 11)

# feature plot
library(ggExtra)
marker_plot=function(SeuratObj,marker){
  pn=length(marker)
  pp=list()
  for(i in 1:pn){
    pg=FeaturePlot(object = SeuratObj, features = marker[i], cols = c("grey", "red"), reduction = "umap")+NoLegend()+labs(x="",y="")+theme(plot.title = element_text(hjust = 0.5))
    pp[[marker[i]]]=pg
  }
  return(pp)
}
p=marker_plot(sce,genes_to_check)
for(i in 1:11){
  ggsave(paste0("plot_out/PNG/all_cell_",genes_to_check[i],"_umap.png"),p[[genes_to_check[i]]],width =6,height =6)
}


#####注释
celltype=data.frame(ClusterID=0:10,
                    celltype= 0:10) 


# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c(0,3,8,9 ),2]='B2'
celltype[celltype$ClusterID %in% c(1,2,4,5,6,10,7),2]='B1'

table(sce@meta.data$RNA_snn_res.0.5)
table(celltype$celltype)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

# rationplot
source('/home/adore_org/Breg/scRNA_scripts/myrationplot.R')
Idents(sce)<-sce$celltype
myrationplot(seurat_object = sce,filename = 'B_type_2',col = mycolors)

## Dimplot_box
source('../scRNA_scripts/Bottom_left_axis.R')
source('../scRNA_scripts/mydimplot.R')
mydimplot(seurat_object = sce,filename = 'B1_B2_umap',reduction = 'umap',group.by = 'celltype',cols = mycolors)

# save
saveRDS(sce, "mymodule_B_celltype.rds")

setwd('../')


source('/home/adore_org/Breg/scRNA_scripts/mymodule_score.R')

B1_marker <- c('Zbtb32', 'Bhlhe41', 'Anxa2', 'Fcrl5', 'Zcwpw1', 'Ass1', 'Zbtb20','Cd9')

B2_marker <- c('Ighd','Fcer2a','Ms4a4c','Vpreb3','Cr2','Ccr7','Neurl3','Icosl')
B1a_marker <- c('Cd19','Cd5','Il10','Itgam')
B1_marker <- marker_class1$B1
B2_marker <- marker_class1$B2
library(stringr)  
B1_marker=str_to_title(B1_marker)
B2_marker=str_to_title(B2_marker)

B1_marker <- B1_marker[(B1_marker %in%  rownames(sce))==TRUE]
B2_marker <- B2_marker[(B2_marker %in%  rownames(sce))==TRUE]
B1a_marker <- B1a_marker[(B1a_marker %in%  rownames(sce))==TRUE]



mymodule_score(RNA,genelist = B1a_marker,cutoff=1.0,filename="B1a_marker")
dir.create('mymodule2')
setwd('mymodule2')
mymodule_score(sce,genelist = B1_marker,cutoff=0.25,filename="B1_marker")
setwd('../')

#####细胞生物学命名

score_B1 <- read.csv('/home/adore_org/B_scRNA-seq/analysis/16-B1/mymodule2/B1_marker_score.csv',header = T,row.names = 1)
score_B2 <- read.csv('/home/adore_org/B_scRNA-seq/analysis/12-B/mymodule2/B2_marker_score.csv',header = T,row.names = 1)
score_B1a <- read.csv('/home/adore_org/B_scRNA-seq/analysis/12-B/mymodule2/B1a_marker_score.csv',header = T,row.names = 1)


#sce <- AddMetaData(sce, metadata = score_B1)
sce$score_cluster_B1 <- score_B1$score_cluster
sce$score_cluster_B2 <- score_B2$score_cluster
sce$score_cluster_B1a <- score_B1a$score_cluster

sce$B1_type <- ifelse(sce$score_cluster_B1=='Positive','B1','NonB1')
sce$B_type <- ifelse(sce$B_type=='B1','B1',ifelse(sce$score_cluster_B2=='Positive','B2','NonB1_B2'))
sce$B1a_type <- ifelse(sce$score_cluster_B1a=='Positive','B1a','NonB1a')
dir.create('B_analysis')
setwd('B_analysis')
#sce$celltype.group <- paste(sce$Breg_type, sce$Sample, sep = "_")
#sce$celltype.group <- factor(sce$celltype.group, levels = c('NonBreg_2MWT','Breg_2MWT','NonBreg_15MWT','Breg_15MWT','NonBreg_26MWT','Breg_26MWT'))
sce$celltype.group <- paste(sce$B_type, sce$Sample, sep = "_")
Idents(sce) <- "celltype.group"

sce$Sample <- factor(sce$Sample, levels=c('26MWT','15MWT','2MWT'))

#
saveRDS(sce,'all_cell_B.rds')

cellfordeg<-levels(sce$Sample)
dir.create('cop_sanmple')
setwd('cop_sanmple')
for(i in 1:length(cellfordeg)){
  for(j in i:length(cellfordeg)){
    CELLDEG <- FindMarkers(sce, ident.1 = paste0('B1_',cellfordeg[i]), ident.2 = paste0('B1_',cellfordeg[j]), only.pos = T,verbose = T,min.pct = 0.25)
    write.csv(CELLDEG,paste0(cellfordeg[i],'_',cellfordeg[j],".CSV"))
}}
setwd('../')
getwd()
list.files()



