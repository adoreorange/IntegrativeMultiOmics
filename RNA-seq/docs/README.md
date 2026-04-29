# RNA-seq 分析流程

## 项目简介

本项目提供了一套完整的 RNA-seq 数据分析流程，包括从原始数据质量控制到差异表达分析和功能富集分析的全过程。

## 分析流程

```
原始数据 (FASTQ)
    ↓
1. 质量控制 (FastQC)
    ↓
2. 数据预处理 (Trimmomatic)
    ↓
3. 序列比对 (STAR)
    ↓
4. 基因定量 (featureCounts)
    ↓
5. 差异表达分析 (DESeq2)
    ↓
6. 功能富集分析 (GO/KEGG/GSEA)
    ↓
7. 结果可视化
```

## 目录结构

```
RNA-seq/
├── code/                    # 分析脚本
│   ├── 01_qc/              # 质量控制
│   │   └── fastqc.sh
│   ├── 02_preprocess/      # 数据预处理
│   │   └── trimmomatic.sh
│   ├── 03_alignment/       # 序列比对
│   │   └── star_align.sh
│   ├── 04_quantification/  # 基因定量
│   │   └── featurecounts.sh
│   ├── 05_de_analysis/     # 差异表达分析
│   │   └── deseq2_analysis.R
│   ├── 06_enrichment/      # 富集分析
│   │   ├── enrichment_analysis.R
│   │   ├── OF_YF_WT_GO.R
│   │   ├── OF_YF_WT_GSEA_analysis.R
│   │   ├── OM_YM_WT_GO.R
│   │   └── TF_enrich.R
│   ├── 07_visualization/   # 可视化
│   │   ├── visualization.R
│   │   ├── GO_KEGG_DOWN.R
│   │   ├── GO_KEGG_UP.R
│   │   ├── GSEA_plot.R
│   │   ├── heatmap_tpm_OF_YF_WT.R
│   │   ├── heatmap_tpm_OM_YM_WT.R
│   │   ├── vacoplot_OF_WT_YF_WT.R
│   │   └── vacoplot_OM_WT_YM_WT.R
│   └── main.sh             # 主流程脚本
├── data/                   # 数据目录
│   ├── raw/               # 原始数据
│   └── sample_info.txt    # 样本信息
├── reference/             # 参考基因组
│   ├── STAR_index/       # STAR 索引
│   └── genes.gtf         # 基因注释文件
├── output/                # 分析结果
├── docs/                  # 文档
└── renv.lock             # R 包版本锁定
```

## 脚本功能说明

### 质量控制 (01_qc/)
- `fastqc.sh`：使用 FastQC 对原始 FASTQ 文件进行质量评估

### 数据预处理 (02_preprocess/)
- `trimmomatic.sh`：使用 Trimmomatic 去除接头序列和低质量碱基

### 序列比对 (03_alignment/)
- `star_align.sh`：使用 STAR 将测序 reads 比对到参考基因组

### 基因定量 (04_quantification/)
- `featurecounts.sh`：使用 featureCounts 对基因进行表达量计数

### 差异表达分析 (05_de_analysis/)
- `deseq2_analysis.R`：使用 DESeq2 进行差异表达分析

### 富集分析 (06_enrichment/)
- `enrichment_analysis.R`：通用富集分析脚本
- `OF_YF_WT_GO.R`：OF vs YF WT 对比组的 GO 富集分析
- `OF_YF_WT_GSEA_analysis.R`：OF vs YF WT 对比组的 GSEA 分析
- `OM_YM_WT_GO.R`：OM vs YM WT 对比组的 GO 富集分析
- `TF_enrich.R`：转录因子富集分析

### 可视化 (07_visualization/)
- `visualization.R`：主可视化脚本
- `GO_KEGG_DOWN.R`：下调基因的 GO/KEGG 富集图
- `GO_KEGG_UP.R`：上调基因的 GO/KEGG 富集图
- `GSEA_plot.R`：GSEA 分析结果可视化
- `heatmap_tpm_OF_YF_WT.R`：OF vs YF WT 对比组的表达量热图
- `heatmap_tpm_OM_YM_WT.R`：OM vs YM WT 对比组的表达量热图
- `vacoplot_OF_WT_YF_WT.R`：OF vs YF WT 对比组的火山图
- `vacoplot_OM_WT_YM_WT.R`：OM vs YM WT 对比组的火山图

## 环境配置

### 依赖软件

| 软件 | 版本 | 用途 |
|------|------|------|
| FastQC | >= 0.11.9 | 质量控制 |
| Trimmomatic | >= 0.39 | 数据预处理 |
| STAR | >= 2.7.10a | 序列比对 |
| featureCounts | >= 2.0.1 | 基因定量 |
| samtools | >= 1.12 | BAM 文件处理 |
| MultiQC | >= 1.11 | QC 报告汇总 |

### R 包依赖

```r
# 差异表达分析
DESeq2
apeglm
tidyverse

# 富集分析
clusterProfiler
org.Hs.eg.db  # 或其他物种注释包
enrichplot
DOSE

# 可视化
ggplot2
pheatmap
RColorBrewer
ggrepel
ComplexHeatmap
patchwork
```

### 安装依赖

#### Conda 环境

```bash
# 创建环境
conda create -n rnaseq python=3.9

# 安装软件
conda install -c bioconda fastqc trimmomatic star subread samtools multiqc

# 安装 R 和 R 包
conda install -c conda-forge r-base r-essentials
```

#### R 包安装

```r
# CRAN 包
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "patchwork"))

# Bioconductor 包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db",
    "enrichplot", "DOSE", "ComplexHeatmap"
))
```

## 使用方法

### 1. 准备数据

将原始 FASTQ 文件放入 `data/raw/` 目录，并创建样本信息文件：

```tsv
# data/sample_info.txt
sample_id    condition    batch
sample1      control      batch1
sample2      control      batch1
sample3      treatment    batch1
sample4      treatment    batch1
```

### 2. 准备参考基因组

```bash
# 下载基因组序列和注释
# 构建 STAR 索引
STAR --runMode genomeGenerate \
     --genomeDir reference/STAR_index \
     --genomeFastaFiles reference/genome.fa \
     --sjdbGTFfile reference/genes.gtf \
     --runThreadN 8
```

### 3. 运行分析

```bash
# 运行完整流程
bash code/main.sh

# 运行指定步骤
bash code/main.sh config/parameters.conf 3 7  # 从步骤3到步骤7

# 运行单个步骤
bash code/01_qc/fastqc.sh data/raw output/01_qc
```

### 4. 配置参数

创建配置文件 `config/parameters.conf`：

```bash
# 数据目录
DATA_DIR="data/raw"
OUTPUT_DIR="output"
REFERENCE_DIR="reference"

# 参考文件
GTF_FILE="${REFERENCE_DIR}/genes.gtf"
STAR_INDEX="${REFERENCE_DIR}/STAR_index"

# 计算资源
THREADS=8
```

## 输出结果

### 1. 质量控制 (01_qc/)
- `fastqc_report.html`: FastQC 质量报告
- `multiqc_report.html`: 汇总 QC 报告

### 2. 比对结果 (03_alignment/)
- `*.sorted.bam`: 排序后的 BAM 文件
- `*_gene_counts.tab`: 基因计数
- `alignment_summary.csv`: 比对统计

### 3. 基因定量 (04_quantification/)
- `count_matrix.txt`: 基因表达矩阵
- `gene_counts.txt`: 详细计数结果

### 4. 差异表达 (05_de_analysis/)
- `DEG_results_all.csv`: 所有基因结果
- `DEG_significant.csv`: 显著差异基因
- `volcano_plot.png`: 火山图
- `PCA_plot.png`: PCA 图
- `heatmap_top50.png`: 热图

### 5. 富集分析 (06_enrichment/)
- `GO_enrichment_*.csv`: GO 富集结果
- `KEGG_enrichment_*.csv`: KEGG 富集结果
- `GSEA_results_*.csv`: GSEA 分析结果
- `TF_enrichment_*.csv`: 转录因子富集结果
- `*_barplot.png`: 柱状图
- `*_dotplot.png`: 气泡图

### 6. 可视化 (07_visualization/)
- `heatmap_tpm_*.png`: TPM 表达量热图
- `volcano_plot_*.png`: 火山图
- `GO_KEGG_*.png`: GO/KEGG 富集图
- `GSEA_plot_*.png`: GSEA 结果图

## 常见问题

### Q: 内存不足怎么办？

A: 调整 STAR 参数或使用 `--genomeLoad LoadAndKeep` 选项共享基因组。

### Q: 如何更改物种？

A: 修改 `enrichment_analysis.R` 中的 `org_db` 和 `kegg_all` 的 organism 参数。

### Q: 如何添加批次效应校正？

A: 在 DESeq2 设计公式中添加 batch 变量：`design = ~ batch + condition`。

## 参考资料

- [DESeq2 文档](https://bioconductor.org/packages/DESeq2)
- [STAR 手册](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [clusterProfiler 教程](https://yulab-smu.top/biomedical-knowledge-mining-book/)

## 作者

RNA-seq 分析流程

## 更新日志

- 2024-01: 初始版本

