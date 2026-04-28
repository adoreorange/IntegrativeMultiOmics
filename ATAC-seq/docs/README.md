# ATAC-seq 分析流程

## 项目简介

本项目提供了一套完整的 ATAC-seq (Assay for Transposase-Accessible Chromatin) 数据分析流程，用于分析染色质开放性数据。

## 分析流程

```
原始数据 (FASTQ)
    ↓
1. 质量控制 (FastQC)
    ↓
2. 数据预处理 (Trimmomatic - Nextera接头去除)
    ↓
3. 序列比对 (BWA)
    ↓
4. Peak Calling (MACS2)
    ↓
5. 差异Peak分析 (DESeq2)
    ↓
6. 功能富集分析 (GO/KEGG)
    ↓
7. 结果可视化
```

## 目录结构

```
ATAC-seq/
├── code/                    # 分析脚本
│   ├── 01_qc/              # 质量控制
│   ├── 02_preprocess/      # 数据预处理
│   ├── 03_alignment/       # 序列比对
│   ├── 04_peak_calling/    # Peak Calling
│   ├── 05_diff_peaks/      # 差异Peak分析
│   ├── 06_enrichment/      # 富集分析
│   ├── 07_visualization/   # 可视化
│   └── main.sh             # 主流程脚本
├── data/                   # 数据目录
│   ├── raw/               # 原始数据
│   └── sample_info.txt    # 样本信息
├── reference/             # 参考基因组
│   ├── BWA_index/        # BWA 索引
│   └── genome.fa         # 基因组序列
├── output/                # 分析结果
├── docs/                  # 文档
└── renv.lock             # R 包版本锁定
```

## 环境配置

### 依赖软件

| 软件 | 版本 | 用途 |
|------|------|------|
| FastQC | >= 0.11.9 | 质量控制 |
| Trimmomatic | >= 0.39 | 数据预处理（Nextera接头去除）|
| BWA | >= 0.7.17 | 序列比对 |
| SAMtools | >= 1.12 | BAM文件处理 |
| MACS2 | >= 2.2.7.1 | Peak Calling |
| bedtools | >= 2.30.0 | BED文件操作 |
| MultiQC | >= 1.11 | QC报告汇总 |

### R 包依赖

```r
# 差异分析
DESeq2
apeglm
tidyverse

# Peak 注释
ChIPseeker
GenomicRanges
rtracklayer

# 富集分析
clusterProfiler
org.Hs.eg.db
TxDb.Hsapiens.UCSC.hg38.knownGene
enrichplot

# 可视化
ggplot2
ComplexHeatmap
pheatmap
RColorBrewer
ggrepel
patchwork
```

### 安装依赖

#### Conda 环境

```bash
# 创建环境
conda create -n atacseq python=3.9

# 安装软件
conda install -c bioconda fastqc trimmomatic bwa samtools macs2 bedtools multiqc

# 安装 R
conda install -c conda-forge r-base r-essentials
```

#### R 包安装

```r
# Bioconductor 包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2", "apeglm", "ChIPseeker", "GenomicRanges",
    "rtracklayer", "clusterProfiler", "org.Hs.eg.db",
    "TxDb.Hsapiens.UCSC.hg38.knownGene", "enrichplot",
    "ComplexHeatmap"
))

# CRAN 包
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "patchwork"))
```

## 使用方法

### 1. 准备数据

将原始 FASTQ 文件放入 `data/raw/` 目录，并创建样本信息文件：

```tsv
# data/sample_info.txt
sample_id    condition    batch
ATAC_rep1    control      batch1
ATAC_rep2    control      batch1
ATAC_rep3    treatment    batch1
ATAC_rep4    treatment    batch1
```

### 2. 准备参考基因组

```bash
# 下载基因组
# 构建 BWA 索引
bwa index -p reference/BWA_index/genome reference/genome.fa
```

### 3. 运行分析

```bash
# 运行完整流程
bash code/main.sh

# 运行指定步骤
bash code/main.sh config/parameters.conf 3 7

# 运行单个步骤
bash code/01_qc/fastqc.sh data/raw output/01_qc
```

## ATAC-seq 特异性处理

### 1. Nextera 接头

ATAC-seq 使用 Tn5 转座酶，接头为 Nextera 接头，而非标准的 Illumina 接头。

### 2. 比对后处理

- **去除重复序列**: 减少 PCR 偏差
- **移除线粒体 reads**: 线粒体 DNA 高度开放，会干扰分析
- **FRiP 分数评估**: Fraction of Reads in Peaks，建议 > 0.3

### 3. Peak Calling 参数

MACS2 ATAC-seq 推荐参数：
```bash
macs2 callpeak \
    -t sample.bam \
    -f BAMPE \
    -g hs \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    --keep-dup all
```

## 输出结果

### 1. 质量控制 (01_qc/)
- FastQC 质量报告
- MultiQC 汇总报告

### 2. 比对结果 (03_alignment/)
- 去重后的 BAM 文件
- 比对统计 (alignment_summary.csv)

### 3. Peak Calling (04_peak_calling/)
- `*_peaks.narrowPeak`: Peak 文件
- `merged_peaks.bed`: 合并的 Peak
- `peak_summary.csv`: Peak 统计
- `frip_scores.csv`: FRiP 分数

### 4. 差异 Peak (05_diff_peaks/)
- `diff_peaks_all.csv`: 所有 Peak 结果
- `diff_peaks_significant.csv`: 显著差异 Peak
- `open_peaks.bed`: 开放 Peak
- `closed_peaks.bed`: 关闭 Peak
- 火山图、MA图、PCA图、热图

### 5. 富集分析 (06_enrichment/)
- Peak 基因组注释
- GO/KEGG 富集结果
- 富集图表

### 6. 可视化 (07_visualization/)
- 染色质开放性热图
- Peak 分布图
- TSS 分布图
- FRiP 分数图

## 质量控制指标

| 指标 | 建议阈值 | 说明 |
|------|----------|------|
| FRiP | > 0.3 | Peak 区域 reads 占比 |
| NFR (Nucleosome Free Region) | 可见 | 片段大小分布应显示 NFR 峰 |
| TSS enrichment | > 10 | TSS 区域富集分数 |
| 线粒体 reads | < 20% | 线粒体 DNA 比例 |

## 常见问题

### Q: 为什么使用 BWA 而不是 Bowtie2？

A: 两者都可以，BWA mem 对 ATAC-seq 的较长片段比对效果较好。

### Q: 如何选择基因组大小参数？

A:
- 人类 (hs): 2.7e9
- 小鼠 (mm): 1.87e9
- 果蝇 (dm): 1.2e8

### Q: FRiP 分数低怎么办？

A:
1. 检查数据质量
2. 调整 peak calling 参数
3. 检查实验质量（细胞状态、转座酶活性）

## 下游分析

### Motif 分析

使用 HOMER 或 MEME 进行 motif 富集分析：

```bash
# HOMER
findMotifsGenome.pl open_peaks.bed hg38 output/motif/

# MEME
# 先提取序列
bedtools getfasta -fi genome.fa -bed open_peaks.bed -fo peaks.fa
meme-chip peaks.fa -o meme_output
```

### 轨道可视化

将结果上传到 IGV 或 UCSC Genome Browser：

```bash
# 转换为 bigWig
bamCoverage -b sample.bam -o sample.bw
bedGraphToBigWig peaks.bedgraph genome.chrom.sizes peaks.bw
```

## 参考资料

- [ENCODE ATAC-seq Pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
- [MACS2 文档](https://macs3-project.github.io/MACS/)
- [ChIPseeker 文档](https://bioconductor.org/packages/ChIPseeker)

## 作者

ATAC-seq 分析流程

## 更新日志

- 2024-01: 初始版本
