# DNA Methylation 分析流程

## 项目简介

本项目提供了一套完整的 DNA 甲基化 (Bisulfite-seq) 数据分析流程，包括从原始数据处理到差异甲基化区域 (DMR) 分析的全过程。

## 分析流程

```
原始数据 (FASTQ)
    ↓
Bismark 比对 (需要单独运行)
    ↓
1. 甲基化提取 (Bismark methylation extractor)
    ↓
2. 甲基化水平计算
    ↓
3. DMR 分析 (DSS)
    ↓
4. 功能富集分析 (GO/KEGG)
    ↓
5. 结果可视化
```

## 目录结构

```
DNA-me/
├── code/                    # 分析脚本
│   ├── 01_preprocess/      # 甲基化提取
│   │   └── bismark_methylation_extractor.sh
│   ├── 02_calculation/     # 甲基化计算
│   │   └── calculate_methylation.R
│   ├── 03_dmr_analysis/    # DMR 分析
│   │   └── dss_dmr_analysis.R
│   ├── 04_enrichment/      # 富集分析
│   │   └── enrichment_analysis.R
│   ├── 05_visualization/   # 可视化
│   │   ├── plot_script/
│   │   │   ├── DMR_female.R
│   │   │   ├── DMR_male.R
│   │   │   ├── YF-WT_vs_OF-WT-ATAC_RRBS-quadrants.R
│   │   │   ├── YF-WT_vs_OF-WT-RRBS_bulk-quadrants.R
│   │   │   ├── YM-WT_vs_OM-WT-ATAC_RRBS-quadrants.R
│   │   │   ├── YM-WT_vs_OM-WT-RRBS_bulk-quadrants.R
│   │   │   ├── mutlAna.R
│   │   │   └── nihe.R
│   │   └── visualization.R
│   └── main.sh             # 主流程脚本
├── data/                   # 数据目录
│   ├── raw/               # 原始数据
│   └── sample_info.txt    # 样本信息
├── reference/             # 参考基因组
│   └── Bismark_genome/   # Bismark 索引
├── output/                # 分析结果
├── docs/                  # 文档
└── renv.lock             # R 包版本锁定
```

## 脚本功能说明

### 甲基化提取 (01_preprocess/)
- `bismark_methylation_extractor.sh`：使用 Bismark 提取甲基化信息

### 甲基化计算 (02_calculation/)
- `calculate_methylation.R`：计算甲基化水平

### DMR 分析 (03_dmr_analysis/)
- `dss_dmr_analysis.R`：使用 DSS 进行差异甲基化区域分析

### 富集分析 (04_enrichment/)
- `enrichment_analysis.R`：DMR 功能富集分析

### 可视化 (05_visualization/)
- `visualization.R`：主可视化脚本
- **plot_script/**：专项分析脚本
  - `DMR_female.R` / `DMR_male.R`：性别特异性 DMR 分析
  - `*quadrants.R`：象限图分析（整合 ATAC-seq 和甲基化数据）
  - `mutlAna.R`：多重分析
  - `nihe.R`：拟合分析

## 环境配置

### 依赖软件

| 软件 | 版本 | 用途 |
|------|------|------|
| Bismark | >= 0.23.0 | 比对和甲基化提取 |
| Bowtie2 | >= 2.4.5 | 比对引擎 |
| SAMtools | >= 1.12 | BAM 文件处理 |

### R 包依赖

```r
# 甲基化分析
methylKit
DSS
bsseq

# 基因组注释
GenomicRanges
rtracklayer
ChIPseeker
annotatr

# 富集分析
clusterProfiler
org.Hs.eg.db
TxDb.Hsapiens.UCSC.hg38.knownGene

# 可视化
ggplot2
ComplexHeatmap
pheatmap
RColorBrewer
patchwork
```

### 安装依赖

#### Conda 环境

```bash
# 创建环境
conda create -n dna-meth python=3.9

# 安装 Bismark
conda install -c bioconda bismark bowtie2 samtools

# 安装 R
conda install -c conda-forge r-base r-essentials
```

#### R 包安装

```r
# Bioconductor 包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "methylKit", "DSS", "bsseq",
    "GenomicRanges", "rtracklayer", "ChIPseeker",
    "clusterProfiler", "org.Hs.eg.db",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "ComplexHeatmap", "annotatr"
))

# CRAN 包
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "patchwork"))
```

## 使用方法

### 1. 准备数据

将原始 FASTQ 文件放入 `data/raw/` 目录，并创建样本信息文件：

```tsv
# data/sample_info.txt
sample_id    condition    batch
BS_rep1      control      batch1
BS_rep2      control      batch1
BS_rep3      treatment    batch1
BS_rep4      treatment    batch1
```

### 2. 准备参考基因组

```bash
# 下载基因组
# 构建 Bismark 索引
bismark_genome_preparation --path_to_bowtie2 /path/to/bowtie2 reference/
```

### 3. 运行 Bismark 比对

```bash
# 单端数据
bismark --genome reference/ -s sample_R1.fq.gz -o output/00_alignment/

# 配对末端数据
bismark --genome reference/ -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
    -o output/00_alignment/ --multicore 8
```

### 4. 运行分析流程

```bash
# 运行完整流程
bash code/main.sh

# 运行指定步骤
bash code/main.sh config/parameters.conf 2 5
```

## 核心概念

### 甲基化水平

DNA 甲基化主要发生在 CpG 二核苷酸的胞嘧啶上。甲基化水平表示为：

```
甲基化水平 (%) = (甲基化 reads / 总 reads) × 100
```

### DMR (Differentially Methylated Region)

差异甲基化区域是两组样本间甲基化水平显著不同的基因组区域。

**分类：**
- **Hyper-methylated**: 处理组甲基化水平显著高于对照组
- **Hypo-methylated**: 处理组甲基化水平显著低于对照组

**筛选标准：**
- FDR < 0.05
- |甲基化差异| > 10%
- 最少 CpG 位点数 ≥ 5

### DML (Differentially Methylated Locus)

差异甲基化位点，指单个 CpG 位点的甲基化差异。

## 输出结果

### 1. 甲基化提取 (01_preprocess/)
- `CX_report.txt`: 甲基化报告
- `*.bedGraph.gz`: 甲基化位点文件
- `methylation_summary.csv`: 甲基化统计

### 2. 甲基化计算 (02_calculation/)
- `methylation_matrix.csv`: 甲基化矩阵
- `sample_methylation_stats.csv`: 样本统计
- `cpg_density_1kb.csv`: CpG 密度

### 3. DMR 分析 (03_dmr_analysis/)
- `DML_all.csv`: 所有差异位点
- `DML_significant.csv`: 显著差异位点
- `DMR_all.csv`: 所有差异区域
- `DMR_significant.csv`: 显著差异区域
- `DMR_hyper.bed` / `DMR_hypo.bed`: BED 文件
- 火山图、MA 图

### 4. 富集分析 (04_enrichment/)
- 基因组注释结果
- GO/KEGG 富集结果
- 富集图表

### 5. 可视化 (05_visualization/)
- 甲基化分布图
- DMR 热图
- PCA 图
- 相关性热图
- 性别特异性 DMR 分析图
- 象限图（整合 ATAC-seq 和甲基化数据）

## 常见问题

### Q: 甲基化水平异常低怎么办？

A: 检查以下几点：
1. BS 转化效率是否正常
2. 文库构建质量
3. 测序深度是否足够

### Q: 如何选择 DMR 筛选阈值？

A: 常用阈值：
- FDR < 0.05（严格）或 0.1（宽松）
- 甲基化差异 > 20% 或 25%
- 根据具体生物学问题调整

### Q: 如何处理批次效应？

A: 在 DSS 分析中加入批次作为协变量：

```r
design <- model.matrix(~ batch + condition)
```

### Q: 如何更改物种？

A: 修改以下参数：
- `org.Hs.eg.db` → `org.Mm.eg.db` (小鼠)
- `TxDb.Hsapiens.UCSC.hg38.knownGene` → `TxDb.Mmusculus.UCSC.mm10.knownGene`
- Bismark 比对时使用相应物种基因组

## 下游分析

### 与基因表达整合

将 DMR 与 RNA-seq 结果整合分析：

```r
# 检查基因启动子区域的甲基化
promoter_dmrs <- annotatePeak(DMR_gr, TxDb = txdb,
                               annoDb = "org.Hs.eg.db",
                               tssRegion = c(-2000, 500))
```

### 与其他表观遗传数据整合

```bash
# 与 ATAC-seq 或 ChIP-seq 数据整合
bedtools intersect -a DMR_hyper.bed -b peaks.bed > overlapping_regions.bed
```

## 参考资料

- [Bismark 文档](https://github.com/FelixKrueger/Bismark)
- [methylKit 文档](https://bioconductor.org/packages/methylKit)
- [DSS 文档](https://bioconductor.org/packages/DSS)
- [DNA Methylation Analysis Guide](https://www.bioconductor.org/packages/release/workflows/vignettes/methylationAnalysisWorkflow/inst/doc/methylationAnalysisWorkflow.html)

## 作者

DNA Methylation 分析流程

## 更新日志

- 2024-01: 初始版本

