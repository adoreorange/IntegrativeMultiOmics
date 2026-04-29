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
dir.create("6-B2")
setwd("6-B2/") 
getwd()

# creat data
set.seed(12345)
sce=readRDS( "/home/adore_org/B_scRNA-seq/analysis/4-B/2-harmony/sce.all_int.rds")
table(sce$RNA_snn_res.1.5)
sce1 = sce[, sce$RNA_snn_res.1.5 %in% c(0,1,3,4,9,10,11,12)]
cellinfo <- subset(sce1@meta.data, select = c("orig.ident", "percent_mito", "Sample", "Barcode",'percent_ribo','percent.redcell'))
sce <- CreateSeuratObject(sce1@assays$RNA@counts, meta.data = cellinfo)

# harmony data

set.seed(10086)
table(sce$orig.ident)
if(T){
  dir.create("2-harmony_B3")
  getwd()
  setwd("2-harmony_B3")
  source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.int = run_harmony(sce)
  setwd('../')
}
sce <- sce.int
sce=readRDS( "./2-harmony/sce.all_int.rds")
sel.clust = "RNA_snn_res.1.5"
sce <- SetIdent(sce, value = sel.clust)

# plot res_tree
p2_tree=clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 9,width = 9)

# plot dimplot
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'umap',group.by = 'RNA_snn_res.1.5');p
ggsave(plot=p, filename="B2_cell_umap_res1.5.pdf",height = 6,width = 7)
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'tsne',group.by = 'RNA_snn_res.1.5');p
ggsave(plot=p, filename="B2_cell_tsne_res1.5.pdf",height = 6,width = 7)
###先根据文中注释看看情况
genes_to_check = c('Zbtb32','Bhlhe41','Fcrl5','Zcwpw1','Cd5', # B1 cells
                   'Ighd','Fcer2a','Ms4a4c','Vpreb3','Cr2','Cd1d1') #B2 ells

library(stringr)  
genes_to_check=str_to_title(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),cols = c("#FFFF00","#FF0000"),
            assay='RNA'  )  + coord_flip();p
ggsave('B2_gene_markers.pdf',height = 7,width = 7)

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
