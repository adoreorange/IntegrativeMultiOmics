#!/usr/bin/env Rscript
rm(list=ls());gc()
setwd('/home/adore_org/B_scRNA-seq/analysis/')
# Run Pseudotime analysis from monocle3

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(monocle3))
library(ggplot2)
library(data.table)
library(stringr)
library(cowplot)
source('/home/adore_org/B_scRNA-seq/analysis/pseudotime/pseudotime_functions.R')
# plot expression umaps
source('/home/adore_org/B_scRNA-seq/analysis/pseudotime/plot_functions.R')



### Parameters ### -------------------------------------------------------------

# specify which data to load
projDir <- '/home/adore_org/B_scRNA-seq/analysis/'
seuDir <- '7-B1/'
seuFile <- './2-harmony_dim15/sce.all_int.rds' # '/Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'
inFile <- paste0(projDir, seuDir, seuFile)

cl_annot <- 'RNA_snn_res.0.7' # which clusters to use
reduction <- 'umap' # 'umap_harmony_t.0.1' # which reduction to use for plotting, and 

rdsDir <- '/home/adore_org/B_scRNA-seq/analysis/9-pseudotime/Rdata/'
plotDir <- '/home/adore_org/B_scRNA-seq/analysis/9-pseudotime/plots/'
output_path <- paste0(rdsDir, seuFile, '_monocle')



### Execute ### ----------------------------------------------------------------
dir.create(plotDir, recursive=T)
dir.create(rdsDir, recursive=T)


# 0. read Seurat file -----
seu = readRDS(inFile)
seu$RNA_snn_res.0.7 <- factor(seu$RNA_snn_res.0.7,levels = c(0,3,8,7,2,9,5,1,4,6))
seu$orig.ident <- factor(seu$orig.ident, levels = c('2MWT','15MWT','26MWT'))
seu$Sample <- factor(seu$Sample, levels = c('2MWT','15MWT','26MWT'))

Idents(seu) <- seu$RNA_snn_res.0.7
m <- seu@meta.data
cl_annot <- 'RNA_snn_res.0.7'
# 6. Run All monocle3 pipeline -----
mon = run.all.monocle3(seu_obj = seu, cl_annot = cl_annot, reduction = reduction)

p <- plot_cells(mon) # plot cells colored by monocle clustering
ggsave(paste0(plotDir,  'trajectory_cluster.pdf'))

p <- plot_cells(mon, color_cells_by='RNA_snn_res.0.7',label_groups_by_cluster = T) # plot cells colored by monocle partition (it'll fit a single graph within each partition)
ggsave(paste0(plotDir, 'trajectory_RNA_snn_res.0.7.pdf'))

p <- plot_cells(mon, color_cells_by='partition') # plot cells colored by cell type name, previously IDd
ggsave(paste0(plotDir, 'trajectory_prediction.pdf'))

p <- plot_cells(mon, color_cells_by='Sample') # plot cells colored by cell type name, previously IDd
ggsave(paste0(plotDir, 'trajectory_Sample.pdf'))

saveRDS(mon, './9-pseudotime/Rdata/monocle.rds')



# # SUBSET the cells used to cells on the pseudotime branch between chol & hep -----
mon <- readRDS(paste0(output_path, '.rds'))

p <- plot_cells(mon, color_cells_by='Sample',
                label_roots=F, label_branch_points = F,
                label_leaves=F, label_cell_groups=F)
p + scale_color_manual(values=c('lightsalmon1', 'skyblue','red'))
p + scale_color_manual(values=mycolors)
p + theme(legend.position='bottom')
ggsave(paste0(plotDir, 'trajectory_cellType.pdf'),
       width=7, height=7)

mon_sub <- choose_graph_segments(mon)
sub_barcodes <- mon_sub@colData@rownames

inSubset <- rownames(colData(mon)) %in% sub_barcodes

# plot the selected cells
mon@colData$selected.Cell <- inSubset

p.selected <- plot_cells(mon, color_cells_by='selected.Cell',
                         label_cell_groups=F,
                         label_branch_points=F,label_roots=F, label_leaves=F)
p.selected.nol <- p.selected+theme(legend.position='none')
p.selected.nol
p.cellType <- plot_cells(mon, color_cells_by='cell.annotation',
                         label_cell_groups=F,
                         label_branch_points=F,label_roots=F, label_leaves=F)
p.cellType
p.cellType.nol <- p.cellType+theme(legend.position='none')

pdf(paste0(plotDir, 'selected_cell_plots.pdf'),
    width=6, height=6)
  print(p.selected)
  print(p.selected.nol)
  print(p.cellType)
  print(p.cellType.nol)
dev.off()



# subset that on the pseudotime segment of interest
seu@meta.data$in.graph.segment <- F
seu@meta.data$in.graph.segment[rownames(seu@meta.data) %in% sub_barcodes] <- T

# subset that also is end stage
keep.barcodes <- rownames(seu@meta.data[seu@meta.data$Disease.status=='end stage' &
                            seu@meta.data$in.graph.segment==T,])

seu@meta.data$keep <- rownames(seu@meta.data) %in% keep.barcodes
Idents(seu) <- seu@meta.data$keep
seu.sub <- subset(seu, idents=TRUE)
seu.sub

# can't work out how to subset monocle a& keep all the stuff needed for plotting
# so need to rerun the monocle pipeline on the subset cells!
mon_sub <- preprocess_cds(mon_sub, method='PCA')
mon_sub <- reduce_dimension(mon_sub, reduction_method='UMAP')
p <- plot_cells(mon_sub)
mon.sub <- run.all.monocle3(seu_obj=seu.sub, cl_annot=cl_annot, reduction=reduction)

p <- plot_cells(mon.sub) # plot cells colored by monocle clustering
ggsave(paste0(plotDir,  'trajectory_cluster_subset.pdf'))

p <- plot_cells(mon.sub, color_cells_by='partition') # plot cells colored by monocle partition (it'll fit a single graph within each partition)
ggsave(paste0(plotDir, 'trajectory_partition_subset.pdf'))

p <- plot_cells(mon.sub, color_cells_by='Disease.status') # plot cells colored by cell type name, previously IDd
ggsave(paste0(plotDir, 'trajectory_Disease_status_subset.pdf'))

p <- plot_cells(mon.sub, color_cells_by='cell.annotation') # plot cells colored by cell type name, previously IDd
ggsave(paste0(plotDir, 'trajectory_cellType_subset.pdf'))


# ORDER CELLS IN PSEUDOTIME ---------
mon.sub <- order_cells(mon.sub)
p <- plot_cells(mon.sub, color_cells_by = 'pseudotime',
                label_leaves=F,
                label_branch_points = F)
p
ggsave(paste0(plotDir, 'trajectory_pseudotime_subset.pdf'))

# Save subset monocle object with trajectory
saveRDS(mon.sub, file = paste0(output_path, '_subset.rds'))
saveRDS(seu.sub, file = paste0(output_path, '_seu_subset.rds'))


# Differential expression analysis ---------------------------------------------
# get genes DE over psuedotime
mon.sub <- readRDS(file=paste0(output_path, '_subset.rds'))

# there is an error in the monocle source code which must be manually edited
# run:
trace('calculateLW', edit=T, where=asNamespace('monocle3'))
# then change "Matrix::rBind" to "rbind". This change isn't permanent, and has to be
# done every time!

# approx 30 mins to run on full
mon.graph.test.res <- graph_test(mon.sub,
                                 neighbor_graph = 'principal_graph',
                                 cores=16)
fwrite(mon.graph.test.res,
       file=paste0(plotDir, 'mon_subset_DE_genes.csv'))
mon.graph.test.res <- fread(paste0(plotDir, 'mon_subset_DE_genes.csv'))



# get genes which have similar expression to 'KFT23', 'FGF13' over the 
# "bridge" ---------------------------------------------------------------------
seu.sub <- readRDS(file=paste0(output_path, '_seu_subset.rds')) # cells in the bridge
mon.sub <- readRDS(file=paste0(output_path, '_subset.rds')) # cells in bridge

# get correlation coeffs for all genes vs the gene of interest
GoI <- 'KRT23'
# # takes a while to run
# corr.df <- get_expression_corr(mon.sub, GoI, corr.method='spearman')
# saveRDS(corr.df, 
#         file=paste0(plotDir, 'KRT23_correlation_df.rds'))
# 
# GoI <- 'FGF13'
# # takes a while to run
# corr.df <- get_expression_corr(mon.sub, GoI, corr.method='spearman')
# saveRDS(corr.df, 
#         file=paste0(plotDir, 'FGF13_correlation_df.rds'))
        

# then plot distributions get the top X genes similar (how many to get to the other goi (FGF13)?) +
# plot them as heatmap (using the plotting functions below)
corr.df.KRT23 <- readRDS(paste0(plotDir, 'KRT23_correlation_df.rds'))
corr.df.FGF13 <- readRDS(paste0(plotDir, 'FGF13_correlation_df.rds'))

# for each of these, make the full table available, and plot the top X genes
corr.df.KRT23 <- corr.df.KRT23[order(corr.df.KRT23$correlation, decreasing=T),]
corr.df.FGF13 <- corr.df.FGF13[order(corr.df.FGF13$correlation, decreasing=T),]

# plot top genes for FGF13 & write full table
# if don't na.omit, get genes with corr = NA, due to 0 variance
DE.genes <- na.omit(corr.df.FGF13$comparison.gene[corr.df.FGF13$correlation > 0.20])
p <- plot_genes_of_interest(mon.sub, DE.genes, font.size=5)
ggsave(paste0(plotDir, 'correlated_genes_FGF13_corr>0.20.pdf'), width=10, height=15)
fwrite(corr.df.KRT23, 
       file=paste0(plotDir, 'correlated_genes_FGF13_allGenes.csv'))

# plot top genes for KRT23 & write full table
# if don't na.omit, get genes with corr = NA, due to 0 variance
DE.genes <- na.omit(corr.df.KRT23$comparison.gene[corr.df.KRT23$correlation > 0.11])
p <- plot_genes_of_interest(mon.sub, DE.genes, font.size=5)
ggsave(paste0(plotDir, 'correlated_genes_KRT23_corr>0.11.pdf'), width=10, height=15)
fwrite(corr.df.KRT23, 
       file=paste0(plotDir, 'correlated_genes_KRT23_allGenes.csv'))




# plotting pseudotime results --------------------------------------------------
seu = readRDS(inFile)
mon.sub <- readRDS(file=paste0(output_path, '_subset.rds'))
seu.sub <- readRDS(file=paste0(output_path, '_seu_subset.rds'))

markers.to.use=3 # choose which cell type markers to plot in c(1,2,3,4)

# 4E like
# extract pseudotime
p <- plot_celltype_pseudotime(seu.sub, mon.sub)
p
ggsave(paste0(plotDir, 'mon_subset_celltype_pseudotime_boxplot.pdf'))


# make analogous plots to fig 4E,F,G,H from Maria Alcolea paper Ilias made
cellType.markers <- fread('/home/USSR/awc30/liver_project/cholangiocyte_hep_analysis/markers/cellType_markers_table_v5.csv')
cellType.markers <- cellType.markers[cellType.markers$marker.gene != 'TFF1']

if (markers.to.use==1) {
  keep.cells <- c('Hepatocytes', 'Cholangiocytes')
  cellType.markers <- cellType.markers[cellType.markers$broad.celltype %in% keep.cells]
} else if (markers.to.use==2) {
  keep.cells <- c('Hepatocytes', 'Cholangiocytes', 
                  'Dedifferentiated progenitor', 'Stem cell')
  cellType.markers <- cellType.markers[cellType.markers$broad.celltype %in% keep.cells]
} else if (markers.to.use==3) {
  keep.cells <- c('Hepatocytes', 'Cholangiocytes', 
                  'Bridge')
  cellType.markers <- cellType.markers[cellType.markers$broad.celltype %in% keep.cells]
} else if (markers.to.use==4) {
  keep.cells <- c('Hepatocytes', 'Cholangiocytes', 
                  'Bridge', 'Dedifferentiated progenitor', 'Stem cell')
  cellType.markers <- cellType.markers[cellType.markers$broad.celltype %in% keep.cells]
} 

# Heatmap of expression by pseudotime with groups for top DE genes
# DE.genes defines the set of genes to plot
DE.genes <- cellType.markers$marker.gene

# get normalised, spline modelled, scaled, subset cells matrix of gene x cells
pt.matrix.sub <- na.omit(prepare_counts_for_plotting(mon.sub,
                                                     DE.genes, 
                                                     N=1))
  
# convert expression matrix to df for plotting
plt.df <- data.frame(pt.matrix.sub)
plt.df$gene <- rownames(plt.df)
plt.dt <- data.table(plt.df)
plt.dt.m <- melt(plt.dt, id.vars='gene')

# cell metadata for plotting
plt.metadata <- make_meta_data_df(mon.sub, seu)

# combine expression, and metadata
plt.dt.m <- merge(plt.dt.m, plt.metadata, by.y='barcode', by.x='variable')

# add marker info about the genes
plt.dt.m <- merge(plt.dt.m, cellType.markers, by.x='gene', by.y='marker.gene', 
                  allow.cartesian=T) # some genes are markers for >1 cellType

# order cells in plot
tmp <- unique(plt.dt.m[, c('variable', 'pseudotime')])
tmp <- tmp[order(tmp$pseudotime),]
plt.dt.m$variable <- factor(plt.dt.m$variable, levels=tmp$variable)

meta.dt.m <- unique(plt.dt.m[, c('variable', 'cell.annotation', 'pseudotime')])
meta.dt.m$variable <- factor(meta.dt.m$variable, levels=tmp$variable)

if (markers.to.use==1) {
  w=7
  h=8
  font.size=20
} else if (markers.to.use %in% c(2, 3, 4)) {
  w=7
  h=8
  font.size=8
}

order.by = 'cellType'
p <- make_pseudotime_expression_heatmap_markers_only(plt.dt.m, meta.dt.m, font.size, order.by=order.by)
p
ggsave(paste0(plotDir,  'mon_subset_pseudotime_heatmap_markersOnly_chol-hep_celltypes_cellTypeOrder_', 
              markers.to.use, '_v5.pdf'), 
       width=w, height=h) 



# Plot UMAPS for the marker genes (those shown in the heatmap) -----------------
gene.vec <- cellType.markers$marker.gene

gv <- gene.vec[1]
pdf(paste0(plotDir, 'slide_14_UMAPs.pdf'), w=10, h=10)
for (gv in sort(gene.vec)) {
  if (gv %in% rownames(seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(seu, gv, is.log=T)
    print(p)
  }
}
dev.off()





