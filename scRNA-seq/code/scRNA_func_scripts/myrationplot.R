library(tidyr)
library(reshape2)
library(ggplot2)
library (gplots) 
library(dplyr)
cols =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
        "#885649","#DD76C5","#BBBE00","#41BED1")
# Barplot of clusters per sample
myrationplot=function(seurat_object,filename,col=cols){
  tb=table((seurat_object@meta.data$orig.ident),(seurat_object@active.ident))
  bar_data <- as.data.frame(tb)
  colnames(bar_data) <- c('Sample', 'Clusters','Freq')
  bar_per <- bar_data %>% group_by(Clusters) %>%
    mutate(sum=sum(Freq)) %>% mutate(percent = Freq / sum)
  tab1 <- cbind(as.data.frame(seurat_object$orig.ident),as.data.frame(seurat_object@active.ident))
  colnames(tab1) <- c("Sample", "Clusters")
  plot1=ggplot(data=tab1,aes(x = Sample, fill = factor(Clusters))) + theme_bw()+theme_cowplot()+
    geom_bar(position = "fill",width=0.6) + theme(axis.text.x = element_text(size = 15,hjust = 1, vjust = 1,angle = 45),axis.text.y = element_text(size = 15),axis.ticks.length = unit(0.3,"cm"))+
    labs(x="",y="Cellular fraction",fill="Cluster")+theme(axis.title.y=element_text(size = 16),legend.title = element_text(size =15)) +
    scale_fill_manual(values=col)
  pdf(paste0(filename,"_all_cell_proportion_per_sample.pdf"),width=12,height=8)
  print(plot1)
  dev.off()
  # Barplot of sample per cluster
  tab2 <- cbind(as.data.frame(seurat_object@active.ident),as.data.frame(seurat_object$orig.ident))
  colnames(tab2) <- c("Clusters", "Sample")
  plot2 = ggplot(data=tab2,aes(x = Clusters, fill = factor(Sample))) +
    geom_bar(position = "fill",width=0.6) + coord_flip()+
    theme(axis.ticks = element_line(linetype = "blank"), legend.position = "top",
          panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
          panel.background = element_rect(fill = NA), plot.background = element_rect(colour = NA)) +
    labs(x="",y="Cellular fraction",fill="Sample")+
    scale_fill_manual(values=col)+ theme_few()+
    theme(plot.title = element_text(size=12,hjust=0.5))+
    #facet_grid(~Var3, scales = "free_x", space ="free_x") +
    theme(plot.title = element_blank(), legend.title = element_blank(),
          legend.key.size = unit(15,"pt"), legend.position = "right",
          axis.text.x = element_text(angle=45, hjust=1, vjust=1))   
  pdf(paste0(filename,"_all_cell_proportion_per_cluster.pdf"),width=10,height=8)
  print(plot2)
  dev.off()
}

