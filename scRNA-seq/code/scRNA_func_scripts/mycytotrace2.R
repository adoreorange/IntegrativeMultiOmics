# 本地安装
# remotes::install_local("./digitalcytometry-cytotrace2-6fe2bad.tar.gz",subdir = "cytotrace2_r", # 特殊的
#                        upgrade = F,dependencies = T)

library(CytoTRACE2)
library(tidyverse)
library(Seurat)

#######输入seurat 对象###########
cytotrace2_result_sce <- cytotrace2(seu, is_seurat = TRUE, slot_type = "counts", 
                                    species = 'mouse',ncores = 15, seed = 1234)
cytotrace2_result_sce

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = seu@meta.data$RNA_snn_res.0.8) %>% set_rownames(., colnames(seu))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, annotation = annotation, is_seurat = TRUE)
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

p5 <- (p1+p2+p3+p4) + plot_layout(ncol = 2)
ggsave(p5, device = "pdf", filename = "./all_CytoTRACE_Plot.pdf", units = "in", width = 12, height = 9)

p6 <- FeaturePlot(MSCs_CT, features = "CytoTRACE_score", reduction = "MSC_UMAP_dim50", 
                  pt.size = 1,coord.fixed = TRUE) + NoAxes() + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p1_raster <- rasterize(p1, dpi = 300)

p6 <- FeaturePlot(cytotrace2_result_sce, reduction = 'umap',features ="CytoTRACE2_Relative",pt.size = 1.5) +
  scale_colour_gradientn(colours = (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)","0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE_score") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + theme(aspect.ratio = 1)
p6_raster <- rasterize(p6, dpi = 300)
ggsave(p6_raster, device = "pdf", filename = "./CytoTRACE_Score_FeaturePlot.pdf", units = "in", width = 5, height = 5)