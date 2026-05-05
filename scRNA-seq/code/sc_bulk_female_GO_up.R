# 绘制sc_bulk_female_GO_up的热图
library(org.Hs.eg.db)
library(clusterProfiler)
source('/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/mycolors.R')
setwd('/home/adore_org/sc_bulk_GO/')
data <- read.csv('./sc_bulk_female_GO.csv',header = T)
use_colors <- c(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                         YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

sc_bulk <- data[which(data$Description %in% c('negative regulation of cell activation','negative regulation of leukocyte activation',
                                              'negative regulation of B cell proliferation','negative regulation of cell population proliferation',
                                              'positive regulation of programmed cell death','negative regulation of cytokine production',
                                              'negative regulation of immune system process','inflammatory response','positive regulation of apoptotic process',
                                              'negative regulation of inflammatory response')),]



# tidy data
sc_bulk <- sc_bulk %>% mutate(log10pvalue=-(sc_bulk$LogP))
sc_bulk <- sc_bulk[order(sc_bulk$log10pvalue, decreasing = T),]
Description <- as.factor(sc_bulk$Description)
sc_bulk$Description <- factor(sc_bulk$Description, levels = rev(Description))

# plot function ----
gk_plot <- ggplot(sc_bulk,aes(Description, y=log10pvalue)) +
  geom_bar(stat="identity", width=0.8, fill='#FF7256') +
  coord_flip() +
  labs(x="", y="- Log10 (P value)", title = 'GO Term Enrichment') +
  theme_pander()  + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line.x = element_line(size = 2, colour = "black"),#x轴连线
        axis.ticks.length.x = unit(0.20, "cm"),#修改x轴刻度的高度，负号表示向上
        axis.text.x = element_text(size = 18, margin = margin(t = 0.3, unit = "cm"), hjust = 0),##线与数字不要重叠 hjust 0，0.5，1
        axis.ticks.x = element_line(colour = "black",size = 1) ,#修改x轴刻度的线    
        axis.title.x = element_text(size = 16), # 修改x轴文本字体大小
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = 16, hjust=1),
        axis.title.y = element_text(size = 16),
        panel.background = element_rect(fill=NULL, colour = 'white'),
        plot.title = element_text(hjust = 0.5, size=16),
        text = element_text(family = "Times") # 调整字体
  ) 
gk_plot

pdf('./enrich_GO_sc_bulk_meta_48.pdf', width = 14, height = 8, useDingbats = T)
gk_plot
dev.off()

