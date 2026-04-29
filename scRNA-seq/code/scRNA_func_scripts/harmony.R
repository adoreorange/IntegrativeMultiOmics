library(Seurat)
library(harmony)
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/Filtervariablefeatures.R')
run_harmony <- function(input_sce){
  print(dim(input_sce))
  input_sce <- NormalizeData(input_sce, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 
  input_sce <- FindVariableFeatures(input_sce,nfeatures=2300)
  input_sce <- Filtervariablefeatures(input_sce, pattern = 'B_cell')
  input_sce <- ScaleData(input_sce)
  print(length(VariableFeatures(object = input_sce)))
  input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
  pdf('ElbowPlot.pdf',width = 8,height = 7)
  pcaplot=ElbowPlot(input_sce, ndims = 40)
  print(pcaplot)
  dev.off()

  seuratObj <- RunHarmony(input_sce, group.by.vars = 'orig.ident', max.iter.harmony = 20,dims.use=1:20)

  names(seuratObj@reductions)
  seuratObj <- RunUMAP(seuratObj,  dims = 1:20, reduction = "harmony")
  #p = DimPlot(seuratObj,reduction = "umap",label=T ) 
  #ggsave(filename='umap-by-orig.ident-after-harmony',plot = p)
  seuratObj <- RunTSNE(seuratObj, dims = 1:20, reduction = "harmony")
  input_sce=seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony", dims = 1:20) 
  input_sce.all=input_sce
  
  #设置不同的分辨率，观察分群效果(选择哪一个？)
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7,0.8, 1, 1.2, 1.5)) {
    input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)
  
  p1_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                     ggtitle("louvain_0.01"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                     ggtitle("louvain_0.1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                     ggtitle("louvain_0.3"))
  
  pdf("Dimplot_diff_resolution_low.pdf",width = 16,height = 8)
  print(p1_dim)
  dev.off()
  p2_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                     ggtitle("louvain_0.8"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                     ggtitle("louvain_1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
                     ggtitle("louvain_0.5"))
  
  pdf("Dimplot_diff_resolution_high.pdf",width = 16,height = 8)
  print(p2_dim)
  dev.off()
  
  p2_tree=clustree(input_sce.all@meta.data, prefix = "RNA_snn_res.")
  ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 9,width = 9)
  table(input_sce.all@active.ident) 
  saveRDS(input_sce.all, "sce.all_int.rds")
  return(input_sce.all)
}