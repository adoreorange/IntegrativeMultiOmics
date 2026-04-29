library(patchwork)
mydimplot=function(seurat_object, filename=filename,reduction = "umap",
                   cols = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282'),
                   pt.size = 0.8, group.by = group.by, label = T,label.box = T){
  pc12 <- Embeddings(object = seurat_object,reduction = reduction) %>%  data.frame()
  # check
  head(pc12,3)
  # get botomn-left coord
  lower <- floor(min(min(pc12$umap_1),min(pc12$umap_2))) - 2
  # get relative line length
  linelen <- abs(0.3*lower) + lower
  # mid point
  mid <- abs(0.3*lower)/2 + lower
  # axies data
  axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                     group = c(1,1,2,2),
                     label = rep(c('umap_2','umap_1'),each = 2))
  # axies label
  label <- data.frame(lab = c('umap_2','umap_1'),angle = c(90,0),
                      x = c(lower - 3,mid),y = c(mid,lower - 2.5))
  # plot
  
  DimPlot<-DimPlot(seurat_object, reduction = reduction, cols = cols, pt.size = 0.8,
                   group.by = group.by,label = T,label.box = T) + NoAxes() + theme(aspect.ratio = 1) + 
    geom_line(data = axes, aes(x = x,y = y,group = group), arrow = arrow(length = unit(0.1, "inches"), ends="last", type="closed")) + 
    geom_text(data = label,aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic') + theme(plot.title = element_blank())+
  pdf(paste0("Dimplot_",filename,"_",reduction, "_box.pdf"),width = 9,height = 7)
  print(DimPlot)
  dev.off()
}
