p <- plot_cells(mon, color_cells_by='RNA_snn_res.0.7',trajectory_graph_color='#476D87',
                label_roots=F, label_branch_points = F,cell_size = 1,group_label_size=6,
                label_leaves=F, label_cell_groups=T) + scale_color_manual(values=mycolors)+ theme(legend.position='right')
pdf('./9-pseudotime/plots/trajectory_RNA_snn_res.0.7.pdf',height = 6,width = 7.5)
p
dev.off()


p <- plot_cells(mon, color_cells_by='Sample',trajectory_graph_color='#476D87',
                label_roots=F, label_branch_points = F,cell_size = 1,group_label_size=4,
                label_leaves=F, label_cell_groups=T) + scale_color_manual(values=mycolors) + theme(legend.position='right')

pdf('./9-pseudotime/plots/trajectory_sample.pdf',height = 6,width = 7)
p
dev.off()

mon <- order_cells(mon)
pdf('./9-pseudotime/plots/pseudotime.pdf',height = 6,width = 7)
plot_cells(mon,
           color_cells_by = "pseudotime",trajectory_graph_color='#476D87',
           label_cell_groups=F,label_roots=F,
           label_leaves=F,cell_size = 1,
           label_branch_points=F,
           graph_label_size=1.5)
dev.off()




#差异基因展示
Track_genes <- graph_test(mon,neighbor_graph="principal_graph", cores=6)
#按莫兰指数选择TOP基因
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
genes <-c('Apoe','Serpinb1a','Cd72','Nfkb1','Tnfrsf13b','Tnfsf8','Ctla4')
plot_genes_in_pseudotime(mon[genes,],
                         color_cells_by="Sample",
                         min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(10)) 

plot_genes_of_interest(mon, genes, font.size=5)

library(monocle)
plot_pseudotime_heatmap(mon[genes,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
