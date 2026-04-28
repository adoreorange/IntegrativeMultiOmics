rm(list=ls());gc()
setwd('/home/adore_org/B_scRNA-seq/analysis/')
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('./scRNA_scripts/mycolors.R')
source('./scRNA_scripts/lib.R')

# creat dir
dir.create("5-B1")
setwd("5-B1/") 
getwd()

# creat data
set.seed(12345)
sce=readRDS("../5-B1/2-harmony/sce.all_int.rds")
table(sce$RNA_snn_res.0.8)
sce1 = sce[, sce$RNA_snn_res.0.8 %in% c(0,1,2,3,4,5,7,8,9)]
cellinfo <- subset(sce1@meta.data, select = c("orig.ident", "percent_mito", "Sample", "Barcode",'percent_ribo','percent.redcell'))
sce1 <- CreateSeuratObject(sce1@assays$RNA@counts, meta.data = cellinfo)

as.data.frame(sce@assays$RNA@counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident) 

# harmony data

set.seed(10086)
table(sce1$orig.ident)
if(T){
  dir.create("2-harmony_dim30")
  getwd()
  setwd("2-harmony_dim30")
  source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.int = run_harmony(sce1)
  setwd('../')
  
}
sce <- sce.int
sce=readRDS("./2-harmony_dim15/sce.all_int.rds")
#sce=readRDS( "./2-harmony/sce.all_int.rds")
sel.clust = "RNA_snn_res.0.7"
sce <- SetIdent(sce, value = sel.clust)

# plot res_tree
p2_tree=clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 9,width = 9)

# plot dimplot
pdf("B_cell_umap_res0.7.pdf",height = 6,width = 7)
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label.size = 5,pt.size = 1,);p
dev.off()


pdf("sample_split_umap.pdf",height = 6,width = 18)
p = DimPlot(sce,label=T,cols = mycolors,reduction = 'umap',split.by = 'orig.ident',label.size = 5,pt.size = 1,);p
dev.off()

p = DimPlot(sce,label=T,cols = mycolors,reduction = 'tsne',group.by = 'RNA_snn_res.0.7');p
ggsave(plot=p, filename="B_cell_tsne_res0.7.pdf",height = 6,width = 7)

source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)
axes <- result$axes
label <- result$label

sample_umap =DimPlot(sce, reduction = "umap",cols = mycolors,pt.size = 0.8,
                     split.by = 'orig.ident',label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label, aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank());sample_umap
ggsave('sample_split_umap.pdf',width = 9,height = 6)


## Dimplot_box
source('../scRNA_scripts/Bottom_left_axis.R')
source('../scRNA_scripts/mydimplot.R')
mydimplot(seurat_object = sce,filename = 'B1',reduction = 'umap',group.by = 'Sample',cols = mycolors)

dir.create('cop_sanmple')
setwd('cop_sanmple')
cellfordeg <- rev(levels(sce$Sample))
for(i in 1:length(levels(sce$Sample))){
  for(j in i:length(levels(sce$Sample))){
    CELLDEG <- FindMarkers(sce, ident.1 = cellfordeg[i], ident.2 = cellfordeg[j], only.pos = F,verbose = T,min.pct = 0.15)
    write.csv(CELLDEG,paste0(cellfordeg[i],'_',cellfordeg[j],".CSV"))
  }}
setwd('../')
getwd()
list.files()

# canlcate marker genes
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/myfindmarkers.R')
dir.create('res0.8')
setwd('res0.8')
Idents(sce) <- sce$RNA_snn_res.0.8
myfindmarkers(sce,filename = 'res0.8',min.pct = 0.1,species = 'mouse',colors = mycolors,thresh.use = 0.25)
setwd('../')
getwd()

#
O_Y <- FindMarkers(sce,ident.1 = c(0,3,5,6),ident.2 = c(1,8),min.pct = 0.15,only.pos = F)

#######
dir.create('cop_sanmple')
setwd('cop_sanmple')
cellfordeg <- rev(levels(sce$Sample))
Idents(sce) <-sce$Sample
for(i in 1:length(levels(sce$Sample))){
  for(j in i:length(levels(sce$Sample))){
    CELLDEG <- FindMarkers(sce, ident.1 = cellfordeg[i], ident.2 = cellfordeg[j], only.pos = F,verbose = T,min.pct = 0.15)
    write.csv(CELLDEG,paste0(cellfordeg[i],'_',cellfordeg[j],".CSV"))
  }}
setwd('../')
getwd()
list.files()

###根据marker注释细胞
genes_to_check = c('Zbtb32','Bhlhe41','Zcwpw1', # B1 cells
                   'Ighd','Fcer2a','Vpreb3')
gene <- c('Apoe','Fos','Jun','Zcwpw1','Cd72', 'Ltb', 'Ly6k', 'Bcl11a','Ctla4')
# features plot
dir.create('featureplot')
setwd('featureplot')
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/myfeatureplot.R')
myfeatureplot(sce,genename = gene,filename = 'Age_gene')
setwd('../')
getwd()

# stacked_violin_plot
dir.create('stack_vlnplot')
setwd('stack_vlnplot')
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/stacked_violin_plot.R')
gene <- c('Apoe','Fos','Fcgr2b','Junb','Jun','Wfdc17','Cyp4f18','Fcrl5','Zcwpw1','Serpinb1a','Rbm3','Nfkb1')
AP <- c('Apoe','Fos','Jun','Junb')
stacked_violin_plot(gene = gene,seurat_object = sce,text.size = 10,flip = F,
                    filename = "Age_gene",width = 12,height = 10,limits.max = 9,Mean = F,col = mycolors)

p1 <- VlnPlot(sce, features = gene2,stack=T,flip = T,
              pt.size=0,  group.by = 'RNA_snn_res.0.7',
              cols=my36colors +theme(legend.position = "none"));p1
ggsave(filename = '28_gene_vlnplot.pdf',plot = p1,height = 18,width = 10)

source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/Vlnplot.R')
Old_gene <-c('Apoe','Fos','Jun','Junb', 'Fcgr2b', 'Fcrl5', 'Ly6k')
Young_gene <-c('Hsph1','Ddit4','Lifr','Sfn','Ckb','Lgals3')
Bcheck_gene <- factor(c('Zbtb32','Bhlhe41','Zcwpw1', 'Fcrl5','Ighd','Fcer2a',
                        'Apoe','Fos','Jun','Junb','Ly6k','Cd72'),levels = c('Zbtb32','Bhlhe41','Zcwpw1', 'Fcrl5','Ighd','Fcer2a',
                                                                            'Apoe','Fos','Jun','Junb','Ly6k','Cd72'))
my_stacked_violin_plot(gene = Bcheck_gene,seurat_object = sce,Clusters = 'RNA_snn_res.0.8',filename = "Bcheck_gene",width = 16,height = 16)

setwd('../')
gene <- c('Apoe','Fos','Fcgr2b','Junb','Jun','Wfdc17','Cyp4f18','Fcrl5','Zcwpw1','Serpinb1a','Rbm3','Nfkb1')
AP <- c('Apoe','Fos','Jun','Junb','Jund')
# doplot
p = DotPlot(sce, features=unique(gene), assay = 'RNA',group.by = 'Sample',) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,size = 10),
        axis.text.y = element_text(size = 12))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expressed") )+ #legend
  scale_color_gradientn(colours = c("blue", "red"));p
#ggsave('check_last_markers.pdf',height = 11,width = 11)
ggsave('RNA_snn_res.0.8_AP_markers.pdf',height = 6,width = 7)
ggsave('age_stage_AP_markers.pdf',height = 8,width = 9)

# rationplot
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/myrationplot.R')
sce$RNA_snn_res.0.7 <- factor(sce$RNA_snn_res.0.7,levels = c(0,3,8,7,2,9,5,1,4,6))
colors <- c('0'='#E64A35','1'='#4DBBD4' ,'2'='#6BD66B','3'='#3C5588','4'='#F29F80' ,'5'='#01A187',
            '6'='#8491B6','7'='#91D0C1','8'='#7F5F48','9'='#AF9E85')
sce$orig.ident <- factor(sce$orig.ident, levels = c('2MWT','15MWT','26MWT'))

Idents(sce) <- sce$RNA_snn_res.0.7
myrationplot(seurat_object = sce,filename = 'B1_res0.7',col = colors)

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('B1_check_markers.pdf',height = 10,width = 7)


#####细胞生物学命名
stage=data.frame(ClusterID=0:8,
                 stage= 0:8) 

# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
stage[stage$ClusterID %in% c(2,4,7),2]='Young'
stage[stage$ClusterID %in% c(1,8),2]='Mid'
stage[stage$ClusterID %in% c(0,3,5,6),2]='Old'

table(sce@meta.data$RNA_snn_res.0.8)
table(stage$stage)

sce@meta.data$stage = "NA"
for(i in 1:nrow(stage)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.8 == stage$ClusterID[i]),'stage'] <- stage$stage[i]}
table(sce@meta.data$stage)

###
sce$age_stage <- paste0(sce$Sample,'_',sce$stage)
sce$age_stage <- factor(sce$age_stage,levels = c('26MWT_Old','15MWT_Old','2MWT_Old','26MWT_Mid','15MWT_Mid','2MWT_Mid',
                                                 '26MWT_Young','15MWT_Young','2MWT_Young'))


#######
dir.create('Age_stage')
setwd('Age_stage')
cellfordeg <- unique(sce$stage)
Idents(sce) <-sce$age_stage

for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(sce, ident.1 = paste0('26MWT','_',cellfordeg[i]), ident.2 = paste0('2MWT','_',cellfordeg[i]), only.pos = F,verbose = T,min.pct = 0.15)
  write.csv(CELLDEG,paste0('26MWT_2MWT','_',cellfordeg[i],".CSV"))
}
setwd('../')
getwd()
list.files()




th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)

result <- left_axes(sce)
axes <- result$axes
label <- result$label
celltype_umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.5,
                       group.by = "celltype",label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
celltype_umap
ggsave('myeloid_umap_by_celltype.pdf',width = 9,height = 7)

saveRDS(sce, "B1_age.rds")
s <- data.frame(gene=sce@assays[["RNA"]]@var.features)
write.csv(s, file='variable_gene.csv')
