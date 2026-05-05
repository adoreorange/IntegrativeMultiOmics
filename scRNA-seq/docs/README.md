# scRNA-seq 分析流程

## 项目简介

本项目提供了一套完整的单细胞 RNA 测序 (scRNA-seq) 数据分析流程，包括从数据质控、批次效应校正到细胞分型、差异表达分析和基因调控网络分析的全过程。

## 分析流程

```
0. 原始数据(FASTQ)
    ↓
1. 数据质控与过滤
    ↓
2. 批次效应校正 (Harmony)
    ↓
3. 细胞聚类与分型
    ↓
4. 细胞亚群分析
    ↓
5. 差异表达基因分析
    ↓
6. 功能富集分析
    ↓
7. 拟时序分析 (Monocle3)
    ↓
8. 基因调控网络分析 (SCENIC)
    ↓
9. 结果可视化
```

## 目录结构

```
scRNA-seq/
├── code/                    # 分析脚本
│   ├── scRNA_func_scripts/ # 功能脚本库
│   │   ├── qc.R
│   │   ├── harmony.R
│   │   ├── myGSVA.R
│   │   ├── myfindmarkers.R
│   │   ├── mydimplot.R
│   │   ├── myfeatureplot.R
│   │   ├── mymodule_score.R
│   │   ├── mycytotrace.R
│   │   ├── mycytotrace2.R
│   │   ├── my_vlnplot_sina.R
│   │   ├── stacked_violin_plot.R
│   │   ├── Heat_dot_data.R
│   │   ├── Filtervariablefeatures.R
│   │   ├── Vlnplot.R
│   │   ├── check-all-markers.R
│   │   ├── lib.R
│   │   ├── myrationplot.R
│   │   ├── Bottom_left_axis.R
│   │   ├── mycolors.R
│   │   └── venn.R
│   ├── 0-cellranger_count.sh    # Cell Ranger 原始比对
│   ├── 1-QC&2-harmony.R    # 质控与批次校正
│   ├── 3-celltype.R        # 细胞类型鉴定
│   ├── 4-B.R               # 细胞亚群分析
│   ├── 5-B1.R              # B1 亚群分析
│   ├── 6-B2.R              # B2 亚群分析
│   ├── 7-DEG.R             # 差异表达分析
│   ├── 8-pyscenic/         # SCENIC 分析
│   │   ├── R_packages.R
│   │   ├── aucell.sh
│   │   ├── csv2loom.py
│   │   ├── get_count_from_seurat.R
│   │   ├── pyscenic_plot.R
│   │   └── run_pyscenic.sh
│   ├── 9-pseudotime/       # 拟时序分析
│   │   ├── monocle3_plot.R
│   │   ├── plot_functions.R
│   │   ├── pseudotime_based_analysis.R
│   │   └── pseudotime_functions.R
│   ├── 10-BCR/             # 免疫组库分析
│   │   ├── B1_immunarch.R
│   │   ├── B1_sample_plot.R
│   │   ├── B2_immunarch.R
│   │   └── B2_sample_plot.R
│   ├── 11-GSVA/            # GSVA 富集分析
│   │   ├── dords.R
│   │   └── gsva.R
│   ├── 12-heatmap/         # 热图可视化
│   │   ├── Heat_Dot_data.R
│   │   └── heatmap_plot.R
│   └── sc_bulk_female_GO_up.R  # 拟 bulk 分析
├── data/                   # 数据目录
│   ├── raw/               # 原始数据
│   └── sample_info.txt    # 样本信息
├── reference/             # 参考数据
│   └── motifs/            # Motif 数据库
├── output/                # 分析结果
├── docs/                  # 文档
└── renv.lock             # R 包版本锁定
```

## 脚本功能说明

### 主分析流程

| 脚本                        | 功能描述                         |
| --------------------------- | -------------------------------- |
| `0-cellranger_count.sh`     | Cell Ranger 原始数据比对         |
| `1-QC&2-harmony.R`          | 数据质控与批次效应校正 (Harmony) |
| `3-celltype.R`              | 细胞类型鉴定                     |
| `4-B.R`, `5-B1.R`, `6-B2.R` | 细胞亚群分析                     |
| `7-DEG.R`                   | 差异表达基因分析                 |
| `sc_bulk_female_GO_up.R`    | 结合 bulk GO 上调分析            |

### SCENIC 分析 (8-pyscenic/)

| 脚本                      | 功能描述                   |
| ------------------------- | -------------------------- |
| `run_pyscenic.sh`         | SCENIC 主分析脚本          |
| `aucell.sh`               | AUCell 分析                |
| `csv2loom.py`             | CSV 转 Loom 格式           |
| `get_count_from_seurat.R` | 从 Seurat 对象提取计数矩阵 |
| `pyscenic_plot.R`         | SCENIC 结果可视化          |
| `R_packages.R`            | R 包加载脚本               |

### 拟时序分析 (9-pseudotime/)

| 脚本                                         | 功能描述        |
| -------------------------------------------- | --------------- |
| `pseudotime_based_analysis.R`                | 拟时序主分析    |
| `monocle3_plot.R`                            | Monocle3 可视化 |
| `pseudotime_functions.R`, `plot_functions.R` | 辅助函数        |

### 免疫组库分析 (10-BCR/)

| 脚本               | 功能描述            |
| ------------------ | ------------------- |
| `B1_immunarch.R`   | B1 细胞免疫组库分析 |
| `B1_sample_plot.R` | B1 样本可视化       |
| `B2_immunarch.R`   | B2 细胞免疫组库分析 |
| `B2_sample_plot.R` | B2 样本可视化       |

### GSVA 分析 (11-GSVA/)

| 脚本      | 功能描述      |
| --------- | ------------- |
| `gsva.R`  | GSVA 富集分析 |
| `dords.R` | 秩相关分析    |

### 热图可视化 (12-heatmap/)

| 脚本              | 功能描述       |
| ----------------- | -------------- |
| `heatmap_plot.R`  | 热图可视化     |
| `Heat_Dot_data.R` | 热图与点图数据 |

### 功能脚本库 (scRNA_func_scripts/)

| 脚本                       | 功能描述         |
| -------------------------- | ---------------- |
| `qc.R`                     | 质量控制         |
| `harmony.R`                | Harmony 批次校正 |
| `myGSVA.R`                 | GSVA 分析        |
| `myfindmarkers.R`          | 差异标记基因查找 |
| `mydimplot.R`              | 降维可视化       |
| `myfeatureplot.R`          | 特征表达图       |
| `mycytotrace.R`            | CytoTRACE 分析   |
| `mycytotrace2.R`           | CytoTRACE 分析2  |
| `my_vlnplot_sina.R`        | Sina 图可视化    |
| `stacked_violin_plot.R`    | 堆叠小提琴图     |
| `mymodule_score.R`         | 模块评分         |
| `myrationplot.R`           | 比例图           |
| `check-all-markers.R`      | 检查所有标记基因 |
| `Vlnplot.R`                | 小提琴图         |
| `Filtervariablefeatures.R` | 筛选可变特征     |
| `Heat_dot_data.R`          | 热图点图数据     |
| `Bottom_left_axis.R`       | 左下坐标轴       |
| `mycolors.R`               | 颜色配置         |
| `lib.R`                    | 库加载           |
| `venn.R`                   | Venn 图          |

## 环境配置

### 依赖软件

| 软件       | 版本     | 用途          |
| ---------- | -------- | ------------- |
| Python     | >= 3.8   | PySCENIC 分析 |
| loompy     | >= 3.0.6 | Loom 文件处理 |
| CellRanger | >= 6.0   | 原始数据比对  |

### R 包依赖

```r
# 单细胞分析
Seurat
monocle3
SingleCellExperiment

# 批次校正
harmony

# 差异分析
DESeq2
edgeR

# 富集分析
clusterProfiler
org.Hs.eg.db
GSVA
fgsea

# 可视化
ggplot2
ComplexHeatmap
pheatmap
RColorBrewer
patchwork

# 免疫组库
immunarch

# 其他
tidyverse
data.table
Matrix
```

### 安装依赖

#### Conda 环境

```bash
# 创建环境
conda create -n scrna python=3.8

# 安装 Python 依赖
pip install pyscenic loompy

# 安装 R
conda install -c conda-forge r-base r-essentials
```

#### R 包安装

```r
# Bioconductor 包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "Seurat", "monocle3", "SingleCellExperiment",
    "harmony", "DESeq2", "edgeR",
    "clusterProfiler", "org.Hs.eg.db", "GSVA", "fgsea",
    "ComplexHeatmap", "immunarch"
))

# CRAN 包
install.packages(c("tidyverse", "data.table", "Matrix",
                   "pheatmap", "RColorBrewer", "patchwork"))
```

## 使用方法

### 1. 准备数据

将原始数据放入 `data/raw/` 目录，可以是：

- FASTQ 文件（需要先进行比对和定量）
- 已处理的 Seurat 对象（RDS 格式）
- 计数矩阵

### 2. 运行分析

```bash
# 运行 Cell Ranger 原始数据比对
bash code/0-cellranger_count.sh

# 运行质控和批次校正
Rscript code/1-QC\&2-harmony.R

# 运行细胞分型
Rscript code/3-celltype.R

# 运行细胞亚群分析
Rscript code/4-B.R
Rscript code/5-B1.R
Rscript code/6-B2.R

# 运行差异表达分析
Rscript code/7-DEG.R

# 运行 GSVA 富集分析
Rscript code/11-GSVA/gsva.R

# 运行热图可视化
Rscript code/12-heatmap/heatmap_plot.R

# 运行 SCENIC 分析
bash code/8-pyscenic/run_pyscenic.sh

# 运行拟时序分析
Rscript code/9-pseudotime/pseudotime_based_analysis.R

# 运行免疫组库分析
Rscript code/10-BCR/B1_immunarch.R
Rscript code/10-BCR/B2_immunarch.R
```

## 核心概念

### 细胞质控指标

| 指标         | 建议阈值 | 说明                   |
| ------------ | -------- | ---------------------- |
| nFeature_RNA | > 200    | 每个细胞检测到的基因数 |
| nCount_RNA   | > 800    | 每个细胞的 UMI 计数    |
| percent.mt   | < 10%    | 线粒体基因比例         |

### 细胞聚类

使用 Seurat 的标准流程：

1. Normalization
2. Feature selection
3. Scaling
4. PCA
5. Harmony batch correction
6. UMAP/tSNE
7. Graph-based clustering

### SCENIC 分析流程

1. **基因共表达网络构建**：使用 GENIE3/GRNBoost2
2. **Regulon 推断**：识别转录因子及其靶基因
3. **AUCell 评分**：计算每个细胞的 regulon 活性

## 输出结果

### 1. 质控结果

- `qc_metrics.csv`: 细胞质量指标
- `filtered_cells.rds`: 过滤后的 Seurat 对象
- QC 统计图

### 2. 细胞分型

- `seurat_object.rds`: 完整的 Seurat 对象（含聚类信息）
- `cell_type_annotations.csv`: 细胞类型注释
- UMAP/tSNE 图

### 3. 差异表达分析

- `DEG_results_*.csv`: 差异表达基因结果
- `marker_genes_*.csv`: 标记基因列表
- 火山图、热图

### 4. 功能富集

- `GO_enrichment_*.csv`: GO 富集结果
- `GSVA_scores.csv`: GSVA 评分矩阵
- 富集分析图

### 5. SCENIC 分析

- `regulons.rds`: Regulon 列表
- `aucell_scores.csv`: AUCell 评分
- Regulon 活性热图

### 6. 拟时序分析

- `pseudotime_results.csv`: 拟时序结果
- `trajectory_plot.png`: 细胞轨迹图

### 7. 免疫组库分析

- `immune_repertoire_*.csv`: 免疫组库信息
- V(D)J 基因使用频率

## 常见问题

### Q: 内存不足怎么办？

A:

1. 使用 `Seurat::RenameCells()` 减少对象大小
2. 对大型数据集使用稀疏矩阵
3. 使用 `gc()` 定期清理内存

### Q: 批次效应校正效果不好怎么办？

A:

1. 尝试不同的批次校正方法（Harmony, Seurat v5 integration, Scanorama）
2. 检查批次是否与生物学变量混淆
3. 考虑移除批次效应明显的样本

### Q: 如何选择聚类分辨率？

A:

- 使用 `FindClusters()` 的 `resolution` 参数
- 通常范围在 0.1-1.0 之间
- 根据生物学知识判断聚类是否合理

### Q: SCENIC 分析需要哪些输入？

A:

- 基因表达矩阵（行=基因，列=细胞）
- 转录因子列表
- Motif 数据库（如 cisTarget）

## 下游分析

### 细胞通讯分析

使用 CellChat 或 NicheNet 分析细胞间通讯：

```r
# CellChat
library(CellChat)
cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")
```

### 多组学整合

与 ATAC-seq 或甲基化数据整合：

```r
# 使用 Signac 整合 scATAC-seq
library(Signac)
combined <- RunATACIntegration(
    RNA = seurat_obj,
    ATAC = atac_obj,
    verbose = TRUE
)
```

## 参考资料

- [Seurat 文档](https://satijalab.org/seurat/)
- [Monocle3 文档](https://cole-trapnell-lab.github.io/monocle3/)
- [SCENIC 文档](https://pyscenic.readthedocs.io/)
- [Harmony 文档](https://github.com/immunogenomics/harmony)
- [CellChat 文档](https://github.com/sqjin/CellChat)

## 作者

scRNA-seq 分析流程

## 更新日志

- 2024-01: 初始版本
