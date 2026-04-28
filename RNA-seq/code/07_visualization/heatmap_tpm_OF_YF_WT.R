library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
library(vegan)
rm(list=ls());gc()
setwd('/home/adore_org/Bulk_data_analysis/')

# data
data <- fread('/home/adore_org/B1_analysis/bulk_genes_counts/rsem.merged.gene_tpm.tsv',header = T, stringsAsFactors = F)
OF <- data.frame(of=rev(grep('OF-WT',colnames(data),value = T)))
YF <- data.frame(yf=rev(grep('YF-WT',colnames(data),value = T)))
deg.all_data <- dplyr::select(data, c(gene_id,YF$yf,OF$of))

features_data <- read.csv('/home/adore_org/Bulk_data_analysis/data_out/OF_YF_WT.csv',header = T, stringsAsFactors = F)
up_features <-features_data$symbol[features_data$change=='up']
down_features <-features_data$symbol[features_data$change=='down']
features <- c(up_features, down_features)
print(length(features))

dif_data <- dplyr::filter(deg.all_data, gene_id %in% features)


dif_data <- as.data.frame(dif_data)
rownames(dif_data)<- dif_data$gene_id
dif_data <- dif_data[,-1]
data.1 <- as.matrix(dif_data)
data.1 <- t(scale(t(data.1)))
#data.1 <- decostand(data.1,"standardize",MARGIN = 1)
#data.1 <- as.data.frame(scale(data))
data.1 <- data.1[match(features, rownames(data.1)), ]

library(ComplexHeatmap)
library(circlize) ## color
# row split

row_split = c()
row_split =factor(row_split,levels = c())
# split heatmap


col_split = c(rep('YF_WT',2),rep('OF_WT',2))
col_split =factor(col_split,levels = c('YF_WT','OF_WT'))

# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue" , "#EEEEEE", "red"), space = "RGB")
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#118ab2", "#fdffb6", "#e63946"))
ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row", annotation_name_side = 'top',annotation_name_gp = gpar(fontsize=14),
                       col = list(Marker = c("MHC-I" = "#e76f51", "TIGIT" = "#ffafcc","NKp46"="#0077b6",'NKG2D'="#00b4d8")))

use_colors <- data.frame(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')
ha_top = HeatmapAnnotation(df = data.frame(Type=col_split),which = 'col', annotation_name_side = 'right',
                           col = list(Type = c('YM_WT' ='#009bff','OM_WT' ='#5558c7',"YF_WT"='#FFA500',"OF_WT" = '#FF4500'))) # 'YF'="#ffafcc",'OF'="#00b4d8"



pdf('./plot_out/heatamp_tpm_OF_YF_WT.pdf',width = 8,height = 12)
pheatmap(data.1,cluster_cols=T,cluster_rows = F,
         show_colnames = T,show_rownames = F,
         heatmap_legend_param = list(title = "Expression",legend_height=unit(3, "cm")),
         color = col_fun,angle_col = '0',labels_col =col_split,)
dev.off()


pdf('./plot_out/heatamp_tpm_OF_YF_WT.pdf',width = 8,height = 12)
Heatmap(data.1, col = col_fun,cluster_columns = T,cluster_rows = F,
        show_column_names = F,show_row_names = F,
        column_names_side = 'bottom',row_names_side = "right",
        column_split = col_split,row_title = NULL,column_title_side = 'bottom',
        row_gap = unit(3, "mm"),column_gap =unit(0, "mm"), row_names_gp = gpar(fontsize=16),
        column_names_rot = 0,width = 2,height = 2,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

write.csv(data.1,file = 'All_heatmap_OF_YF_WT.csv')

Heatmap(data.1, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = F,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = colome_split,row_title = NULL,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), row_names_gp = gpar(fontsize=16),
        top_annotation = ha_top, left_annotation = ha,column_names_rot = 0,width = 2,height = 2,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
