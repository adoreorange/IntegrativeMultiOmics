library(ggplot2)
rm(list=ls());gc()
use_colors <- c(YM ='#009bff', OM ='#5558c7',YF ='#FFA500',OF ='#FF4500',
                YM_KO ='#8A2BE2', OM_KO ='#130780', YF_KO ='#FF7256', OF_KO ='#bb0a1e')

use_colors <- c('YM_Opened' ='#009bff', 'OM_Opened' ='#5558c7','YF_Opened' ='#FFA500','OF_Opened' ='#FF4500')

peaks <- data.frame(Group=c('OF_Opened','YF_Opened','OM_Opened','YM_Opened'),
                    peaks=c(41176,1352,6900,1116))
peaks$Group <-factor(peaks$Group,levels =c('OF_Opened','YF_Opened','OM_Opened','YM_Opened') )
peaks <- peaks %>% dplyr::mutate(range=peaks/100)

p <- ggplot(data = peaks, aes(x = Group, y = (range), fill = Group)) +
  geom_col() +  # 使用geom_col()来绘制条形图
  scale_fill_manual(values = use_colors) +  # 使用颜色方案
  labs(x='DARs',y='Number of Regions(x100)',title = "") +  # 添加轴标签和图例标题
  theme(legend.title = element_text(size = 14, face = "bold", colour = "black",hjust = 0),
        legend.text = element_text(size = 12, family = 'Times',hjust = 0),
        panel.grid.major = element_blank(),    # 移除主要网格线
        panel.grid.minor = element_blank(),    # 移除次要网格线
        #axis.ticks.x = element_blank(),       # 用于移除x轴的刻度线
        axis.line.x = element_line(size = 1, colour = "black"), # 设置x轴的线条宽度为2，颜色为黑色
        axis.ticks.length.x = unit(0.20, "cm"), # 设置x轴刻度线的长度为0.20cm
        axis.text.x = element_text(size = 14, colour = 'black',margin = margin(t = 0.3, unit = "cm"), hjust = 0.5), # 设置x轴刻度标签的字体大小为18，上边距为0.3cm，水平对齐方式为0左，0.5中，1右
        axis.ticks.x = element_line(colour = "black",size = 1) , # 设置x轴刻度线的颜色为黑色，大小为1    
        axis.title.x = element_text(size = 16), # 设置x轴标题的字体大小为16
       # 用于移除y轴的刻度线
        axis.ticks.y = element_line(colour = "black",size = 1),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text.y  = element_text(size = 18, hjust=1,colour = 'black'), # 设置y轴刻度标签的字体大小为16，水平对齐方式为1（右对齐）
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小为16
        panel.background = element_rect(fill='white', colour = NULL ), # 设置图表面板背景为白色，fill=NULL表示背景不填充颜色
        plot.title = element_text(hjust = 0.5, size=16), # 设置图表主标题的水平对齐方式为居中（0.5），字体大小为16
        text = element_text(family = "Times",colour = 'black'))  # 设置所有文本字体
pdf('./plot/DARs_barplot.pdf',height = 7,width = 7)
p
dev.off()
