library(tidyverse)

# 读取数据
setwd('/home/hyf/rrbs')
my_colors <- c(
  "Promoter" = "#E64B35",
  "5' UTR" = "#4DBBD5",
  "3' UTR" = "#00A087",
  "Exon" = "#3C5488",
  "Intron" = "#F39B7F",
  "Downstream" = "#8491B4",
  "Distal Intergenic" = "#91D1C2",
  "Other" = "#CCCCCC"
)


my_colors <- c(
  "Promoter" = "#D55E00",
  "5'UTR" = "#0072B2",
  "3'UTR" = "#009E73",
  "Exon" = "#CC79A7",
  "Intron" = "#F0E442",
  "Downstream" = "#56B4E9",
  "Intergenic" = "#999999",
  "Other" = "#DDDDDD"
)
dmr <- read.csv("./data/peakAnno_male_dmr.csv", stringsAsFactors = FALSE)


dmr$annotation_clean <- sub(" \\(.*\\)", "", dmr$annotation)
# 标准化 annotation 分类
dmr$annotation_simple <- case_when(
  grepl("Promoter", dmr$annotation_clean) ~ "Promoter",
  grepl("5UTR", dmr$annotation_clean) ~ "5'UTR",
  grepl("3UTR", dmr$annotation_clean) ~ "3'UTR",
  grepl("Exon", dmr$annotation_clean) ~ "Exon",
  grepl("Intron", dmr$annotation_clean) ~ "Intron",
  grepl("Downstream", dmr$annotation_clean) ~ "Downstream",
  grepl("Intergenic", dmr$annotation_clean) ~ "Intergenic",
  TRUE ~ "Other"
)

dmr_stat <- dmr %>%
  group_by(annotation_clean) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# 计算比例
dmr_stat <- dmr_stat %>%
  mutate(percent = count / sum(count) * 100)

dmr_stat

setwd('/home/hyf/rrbs/plot')
library(ggplot2)
##柱状图
ggplot(dmr_stat, aes(x = reorder(annotation_clean, -count), 
                     y = count, 
                     fill = annotation_clean)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  labs(x = "Genomic annotation", 
       y = "Number of DMRs",
       title = "Genomic distribution of DMRs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("./male/DMR_barplot.pdf", width = 6, height = 5)

##饼图
ggplot(dmr_stat, aes(x = "", y = percent, fill = annotation_clean)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = my_colors) +
  theme_void() +
  labs(title = "DMR genomic feature distribution")+
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )
ggsave("./male/DMR_pieplot.pdf", width = 5, height = 5)

dmr_chr <- dmr %>%
  group_by(seqnames) %>%
  summarise(count = n())

chrom_order <- paste0("chr", c(as.character(1:19), "X", "Y"))
dmr_chr$seqnames <- factor( dmr_chr$seqnames, levels = chrom_order)

### 染色体分布
ggplot(dmr_chr, aes(x = seqnames, y = count)) +
  geom_bar(stat = "identity", fill = "#4DBBD5") +
  theme_bw() +
  labs(x = "Chromosome", y = "DMR count",
       title = "DMR distribution across chromosomes")

ggsave("./male/DMR_chr_distribution.pdf", width = 10, height = 5)

## TSS距离（高级灰）
ggplot(dmr, aes(x = distanceToTSS)) +
  geom_histogram(bins = 60, fill = "#3C5488") +
  theme_bw() +
  labs(title = "Distance of DMRs to TSS",
       x = "Distance to TSS", y = "Count")

dmr$region_simple <- ifelse(grepl("Promoter", dmr$annotation_clean),
                            "Promoter", "Non-Promoter")

table(dmr$region_simple)

### Promoter vs Non-promoter（对比色）
ggplot(dmr, aes(x = region_simple, fill = region_simple)) +
  geom_bar() +
  scale_fill_manual(values = c("Promoter" = "#E64B35",
                               "Non-Promoter" = "#4DBBD5")) +
  theme_bw() +
  labs(title = "Promoter vs Non-promoter DMRs")

ggplot(dmr_stat, aes(x = "", y = percent, fill = annotation_clean)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "DMR genomic feature distribution")
