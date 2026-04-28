library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
###### Heatmap plot
setwd('/home/adore_org/B_scRNA-seq/analysis/')
sce <- readRDS('./7-B1/2-harmony_dim15/sce.all_int.rds')
source('../scRNA_scripts/Heat_dot_data.R')
aged <- c('Apoe','Fos','Gapdh','Tnfrsf13c','Nfkb1','Cd9','Hmgb2','Bax','Hnrnpd','Rbm38','Mif','Zfp36','Mdh1','Pkm','Psmd14','Asph','Ewsr1','Ypel3')
Up <- c('Apoe','Serpinb1a','Rbm3','C130026I21Rik','Rcn3','Ms4a6c','Tnfaip8','Cd72','D10Wsu102e',
        'Vars','Cd24a','Tg','AW112010','AI427809','Laptm5','Marcksl1','St3gal6','Nacc2','Zcchc18',
        'Id3','Asph','Nfkb1','Tnfrsf13b','Man1a','Unc119','Gimap3','Cyp11a1','Myo1e','Ssbp3',
        'Scn4a','Gimap7','As3mt','A530032D15Rik','Il5ra','Krt222','Cpt1a','Gimap4','Gm8369',
        'Cdc25b','Hmces','Stk3','Ctla4','Tnfsf8','Ctse','Inpp5k','Cybb','Themis2','Slfn1')
Down <- c('Anxa2','Lifr','Cnn3','Adm','Emb','Pik3ip1',
          'Apobec1','Dnaja1','Emp3','Pdia3','Lgals1','Lgals3','Ddit4',
          'Nrros','Hsp90b1','Tubb2a','Tyrobp','B4galnt1','Plaur','Tubb6','Fgd2','Fam129c',
          'Sfn','Gas7','Hsph1','Chordc1','Sdf2l1','Tmem123','Klf13',
          'Atp8b2','Cxcr4','Akr7a5','Fcer1g','Tubb2b','1700030K09Rik','Zfp36l2')
features = c('Tnfaip8','Serpinb1a','Apoe','Rbm3','Cd72', 'Cd24a','Nacc2','Zcchc18','Id3',
             'Asph','Gimap7','Cdc25b','Slfn1','Nfkb1')
gene <- c('Serpinb1a','Tnfaip8','Cd72','Cd24a','Ctla4','Tnfsf8','Cybb','Il5ra','Gimap3','Gimap4','Gimap7',
          'Tg','Cpt1a','Cyp11a1','Vars','As3mt','Apoe',
          'Nacc2','Cdc25b','Stk3','Hmces',
          'Rbm3','Id3','Nfkb1',
          'Themis2','Inpp5k','Scn4a','Ctse','Myo1e','Ssbp3','Unc119','Marcksl1',
          'Rcn3','Asph','St3gal6','Ms4a6c','Laptm5','Krt222','Man1a')
features = Up
row_split = c(rep('immune-inflammatory response',11),rep('Metabolism and Hormones',6),rep('Cell growth and differentiation',4),
              rep('Transcription and regulation',3),rep('signal transduction',4),rep('Cell Structure and Transport',4),rep('Glycosylation and protein processing',7))
row_split =factor(row_split,levels = c('immune-inflammatory response','Metabolism and Hormones','Cell growth and differentiation',
                                       'Transcription and regulation','signal transduction','Cell Structure and Transport','Glycosylation and protein processing'))

### plot heatmap
source("../scRNA_scripts/Heat_dot_data.R")
### set colnames order
plot_ord <- levels(sce$RNA_snn_res.0.7)
data.plot <- Heat_Dot_data(object=sce,features=features,group.by="RNA_snn_res.0.7")
exp.mat <- data.plot %>% dplyr::select(features.plot,id,avg.exp.scaled) %>% tidyr::spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% dplyr::select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100
write.csv(per.mat,'./heatdata/sc_bulk_up_39_anno.csv')

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color
# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue" , "#EEEEEE", "red"), space = "RGB")
col_fun <- colorRamp2(c(-1.5, 0, 2), c("#118ab2", "#fdffb6", "#e63946"))
# split heatmap
col_split = c(rep('Young',3),rep('Mid',3),rep('Old',4))
col_split =factor(col_split,levels =  c('Young','Mid',"Old"))
# left annotation
annot = c('immune-inflammatory response','Metabolism and Hormones','Cell growth and differentiation','Transcription and regulation','signal transduction','Cell Structure and Transport','Glycosylation and protein processing')
row_color=c("#e76f51","#ffafcc","#0077b6","#ddbea9","#00b4d8","#dc2f02","#2a9d8f","#57cc99",'#7294D4')

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("immune-inflammatory response" = "#e76f51", "Metabolism and Hormones" = "#ffafcc","Cell growth and differentiation"="#0077b6",
                                             'Transcription and regulation'='#ddbea9','signal transduction'='#00b4d8','Cell Structure and Transport'='#dc2f02',
                                             'Glycosylation and protein processing'='#2a9d8f')))
pdf(paste0("./Heatplot/48_gene_heatmap_anno.pdf"),width = 6,height = 16)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"),
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0("./Heatplot/39_gene_dotplot_anno.pdf"),width = 8,height = 16)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/3 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        column_split = col_split,left_annotation = ha,row_title = F,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"),
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()



Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        column_split = col_split,row_split = row_split,left_annotation = ha,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), row_title = F,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))