#!/bin/bash
# ==========================================================
# 10x单细胞数据Cell Ranger Count处理脚本
# 功能：将10x fastq数据转换为单细胞矩阵（.h5或.mtx格式）
# 依赖：Cell Ranger（需提前安装并配置环境变量）
# ==========================================================

# ------------------------------
# 1. 自定义参数（需根据实际数据修改）
# ------------------------------
REF_GENOME="/path/to/refdata-gex-mm10-2020-A"  # 参考基因组路径（如10x官网下载的mm10）
FASTQ_DIR="/path/to/your/fastq_files"           # 包含fastq文件的目录（需包含Illumina BCL转化的fastq）
OUTPUT_DIR="/path/to/output"                   # 输出目录（Cell Ranger会在此生成结果文件夹）
SAMPLE_NAME="your_sample_id"                   # 样本名（与fastq文件中的样本ID一致）
EXPECTED_CELLS=5000                           # 预期细胞数（根据实验设计调整，如5000-10000）
CPU_CORES=16                                  # 使用的CPU核心数（建议不超过服务器总核心数）
MEMORY_GB=64                                  # 使用的内存（GB，建议≥32GB）

# ------------------------------
# 2. 执行Cell Ranger Count
# ------------------------------
cellranger count \
  --id=${SAMPLE_NAME} \                       # 输出目录名（Cell Ranger会在OUTPUT_DIR下生成该文件夹）
  --transcriptome=${REF_GENOME} \             # 参考基因组路径
  --fastqs=${FASTQ_DIR} \                     # fastq文件所在目录
  --sample=${SAMPLE_NAME} \                   # 样本名（需与fastq文件中的样本ID匹配）
  --expect-cells=${EXPECTED_CELLS} \           # 预期细胞数（提高细胞检测准确性）
  --localcores=${CPU_CORES} \                 # 并行使用的CPU核心数
  --localmem=${MEMORY_GB}                     # 分配的内存（GB）

# ------------------------------
# 3. 脚本说明
# ------------------------------
# 输出结果包含：
# - matrices/h5：单细胞表达矩阵（.h5格式，兼容Seurat/Scanpy）
# - filtered_feature_bc_matrix：过滤后的表达矩阵（.mtx格式，用于下游分析）
# - web_summary.html：分析报告（可视化质控结果）
#
# 注意事项：
# - 确保fastq目录结构符合10x标准（如含`Sample_S1_L001_R1_001.fastq.gz`等文件）
# - 参考基因组需与实验物种/版本匹配（如小鼠用mm10，人类用GRCh38）
# - 若fastq文件跨多个lane，Cell Ranger会自动合并
