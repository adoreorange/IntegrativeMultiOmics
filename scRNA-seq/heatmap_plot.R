library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
###### Heatmap plot
features=c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4",
           "RUNX1","PTPRC","SPN",
           "HLF","GFI1","MYB","MYCN",
           "C1QA","CD14","LYVE1","LYZ","RNASE2",
           "COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1",
           "EPCAM","AFP","FGB","APOA1","MAL","CALB1",
           "HBE1","HBZ","GYPA",
           "MKI67","TOP2A","AURKB")
row_split = c(rep("Endo",3),rep("VEC",2),rep("AEC",2),rep("Hema",3),rep("HSC",4),rep('Mo/Mø',3),rep('Gr',2),rep('Stroma',5),rep('Fibroblasts',11),rep('Peric',2),rep('SM',2),rep('Podoc',3),rep('Ep',2),rep('Liver',2),rep('Kidn',2),rep('Ery',3),rep('Prolif',3))
row_split =factor(row_split,levels = c("Endo","VEC","AEC","Hema","HSC",'Mo/Mø','Gr','Stroma','Fibroblasts','Peric','SM','Podoc','Ep','Liver','Kidn','Ery','Prolif'))

### plot heatmap
source("./Rcode/Heat_Dot_data.R")
### set colnames order
plot_ord <- levels(sample$seurat_clusters)
data.plot <- Heat_Dot_data(object=sample,features=features,group.by="seurat_clusters")
exp.mat <- data.plot %>% dplyr::select(features.plot,id,avg.exp.scaled) %>% tidyr::spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% dplyr::select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color
# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("#118ab2", "#fdffb6", "#e63946"))
# split heatmap
col_split = c(rep('Endo',4),rep('Hema',2),'HSC',rep('Stro/Me',10),rep('Epith',2),'Ery')
col_split =factor(col_split,levels = c("Endo", "Hema", 'HSC', "Stro/Me", "Epith", "Ery"))
# left annotation
annot = c("Endo", "VEC", "AEC", "Hema", "HSC",  "Mo/Mø", "Gr","Stroma", "Fibroblasts", "Peric", "SM","Podoc", "Ep", "Liver","Kidn","Ery" , "Prolif")
row_color=c("#e76f51","#ffafcc","#0077b6","#ddbea9","#00b4d8","#dc2f02","#2a9d8f","#57cc99","#b5838d","#8a5a44","#023047","#75A8BF", "#C54575", "#403044","#C78D42", "#C992B1", "#804133")

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("Endo" = "#e76f51", "VEC" = "#ffafcc","AEC"="#0077b6","Hema"="#ddbea9",
                                             "HSC"="#00b4d8","Mo/Mø"="#dc2f02","Gr"="#fca311","Stroma"="#57cc99","Fibroblasts"="#b5838d","Peric"="#8a5a44","SM"="#023047",'Podoc'="#75A8BF", 'Ep'="#C54575", 'Liver'="#403044",'Kidn'="#C78D42", 'Ery'="#C992B1", 'Prolif'="#804133")))
pdf(paste0("./plot_out/Cell_identity_AGM_heat_new.pdf"),width = 9,height = 13)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), 
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0("./plot_out/Cell_identity_AGM_heat_new2.pdf"),width = 7,height = 15)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()
