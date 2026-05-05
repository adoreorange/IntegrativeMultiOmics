# 多组学整合分析项目

本项目是一个综合性的多组学数据分析框架，涵盖基因组学、转录组学、表观组学等多个层面的数据处理与分析流程。

## 项目结构

```
AnaCode/
├── ATAC-seq/        # 染色质可及性测序分析
├── DNA-me/          # DNA甲基化分析
├── RNA-seq/         # RNA测序分析
├── scRNA-seq/       # 单细胞RNA测序分析
└── README.md        # 项目说明文档
```

## 模块功能说明

### 1. ATAC-seq 分析模块

**功能定位**：分析染色质开放区域，识别转录调控元件

**分析流程**：
| 步骤 | 目录 | 脚本 | 功能描述 |
|------|------|------|----------|
| 01 | qc | `fastqc.sh` | 原始数据质量控制 (FastQC) |
| 02 | preprocess | `trimmomatic.sh` | 序列预处理 (Trimmomatic) |
| 03 | alignment | `bwa_align.sh` | 序列比对 (BWA) |
| 04 | peak_calling | `macs2_peak.sh` | 峰检测 (MACS2) |
| 05 | diff_peaks | `deseq2_peaks.R` | 差异峰分析 (DESeq2) |
| 06 | enrichment | `enrichment_peaks.R`, `motif_analysis.sh` | 富集分析与基序分析 |
| 07 | visualization | `visualization_peaks.R`, `footprint.sh`, `tss_plot.sh` | 结果可视化、足迹分析、TSS图 |
| - | diffbind | `all_diffbind.R`, `ChIPseeker_all_wt.R`, `enrichGO_*.R`, `enrichKEGG_*.R` | DiffBind差异结合分析、GO/KEGG富集 |

**核心脚本**：`code/main.sh` - 一键执行完整分析流程

---

### 2. DNA-me 分析模块

**功能定位**：分析DNA甲基化水平，探索表观遗传调控

**分析流程**：
| 步骤 | 目录 | 脚本 | 功能描述 |
|------|------|------|----------|
| 01 | preprocess | `bismark_methylation_extractor.sh` | 甲基化提取 (Bismark) |
| 02 | calculation | `calculate_methylation.R` | 甲基化水平计算 |
| 03 | dmr_analysis | `dss_dmr_analysis.R` | 差异甲基化区域分析 (DSS) |
| 04 | enrichment | `enrichment_analysis.R` | 富集分析 |
| 05 | visualization | `visualization.R` | 主可视化脚本 |
| - | plot_script | `DMR_female.R`, `DMR_male.R`, `*quadrants.R`, `mutlAna.R`, `nihe.R` | 性别特异性DMR分析、象限图、多重分析 |

**核心脚本**：`code/main.sh` - 一键执行完整分析流程

---

### 3. RNA-seq 分析模块

**功能定位**：分析基因表达水平，揭示转录组动态变化

**分析流程**：
| 步骤 | 目录 | 脚本 | 功能描述 |
|------|------|------|----------|
| 01 | qc | `fastqc.sh` | 原始数据质量控制 (FastQC) |
| 02 | preprocess | `trimmomatic.sh` | 序列预处理 (Trimmomatic) |
| 03 | alignment | `star_align.sh` | 序列比对 (STAR) |
| 04 | quantification | `featurecounts.sh` | 基因表达定量 (featureCounts) |
| 05 | de_analysis | `deseq2_analysis.R` | 差异表达分析 (DESeq2) |
| 06 | enrichment | `enrichment_analysis.R`, `OF_YF_WT_GO.R`, `OF_YF_WT_GSEA_analysis.R`, `OM_YM_WT_GO.R`, `TF_enrich.R` | GO/KEGG/GSEA富集分析、转录因子富集 |
| 07 | visualization | `visualization.R`, `GO_KEGG_DOWN.R`, `GO_KEGG_UP.R`, `GSEA_plot.R`, `heatmap_tpm_*.R`, `vacoplot_*.R` | 火山图、热图、GSEA图 |

**核心脚本**：`code/main.sh` - 一键执行完整分析流程

---

### 4. scRNA-seq 分析模块

**功能定位**：单细胞水平基因表达分析，解析细胞异质性

**分析流程**：
| 步骤 | 脚本 | 功能描述 |
|------|------|----------|
| 0 | `0-cellranger_count.sh` | Cell Ranger 原始输出 |
| 1 | `1-QC&2-harmony.R` | 数据质控与批次效应校正 (Harmony) |
| 3 | `3-celltype.R` | 细胞类型鉴定 |
| 4-6 | `4-B.R`, `5-B1.R`, `6-B2.R` | 细胞亚群分析 |
| 7 | `7-DEG.R` | 差异表达基因分析 |
| 10 | `10-BCR/` | 免疫组库分析 (B1/B2_immunarch.R, B1/B2_sample_plot.R) |
| 11 | `11-GSVA/` | GSVA功能富集分析 |
| 12 | `12-heatmap/` | 热图与点图可视化 |

**子目录说明**：
| 目录 | 内容 |
|------|------|
| `8-pyscenic/` | SCENIC基因调控网络分析 |
| `9-pseudotime/` | 拟时序分析 (Monocle3) |
| `10-BCR/` | B细胞受体免疫组库分析 |
| `11-GSVA/` | GSVA功能富集分析 |
| `12-heatmap/` | 热图可视化 |
| `scRNA_func_scripts/` | 功能脚本库 (qc.R, harmony.R, myGSVA.R, myfindmarkers.R, mydimplot.R, mycytotrace.R等) |

---

## 技术栈

| 类别 | 工具/包 |
|------|---------|
| **测序质控** | FastQC, Trimmomatic |
| **序列比对** | BWA, STAR |
| **表观分析** | MACS2, Bismark |
| **表达定量** | featureCounts |
| **差异分析** | DESeq2, DSS |
| **富集分析** | clusterProfiler, GSVA, GSEA |
| **可视化** | ggplot2, pheatmap, ComplexHeatmap |
| **单细胞分析** | Seurat, Monocle3, SCENIC, immunarch |
| **批次校正** | Harmony |

## 使用方法

每个分析模块都提供了 `main.sh` 脚本，可通过以下方式执行：

```bash
# 进入对应模块目录
cd ATAC-seq/code/

# 执行主脚本
bash main.sh
```

或单独执行各步骤脚本：

```bash
# 以RNA-seq为例，执行质量控制
bash 01_qc/fastqc.sh
```

---

## 目录说明

| 目录/文件 | 用途 |
|-----------|------|
| `code/` | 分析脚本存放目录 |
| `docs/` | 模块详细说明文档 |
| `renv.lock` | R包环境依赖锁定文件 |
| `code/main.sh` | 一键执行完整分析流程 |

---

## 注意事项

1. 运行前请确保已安装所有必需的生物信息学工具和R包
2. 各模块的 `renv.lock` 可用于复现R环境
3. 建议在Linux/Unix环境下运行分析脚本
4. 大型数据集分析需要充足的计算资源和内存

---

**项目维护者**：adoreorange  
**最后更新**：2026年4月

