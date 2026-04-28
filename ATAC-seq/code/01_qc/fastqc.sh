#!/bin/bash
# FastQC 质量控制分析 (ATAC-seq)
# 用法: bash fastqc.sh <input_dir> <output_dir>

set -e

# 参数设置
INPUT_DIR=${1:-"data/raw"}
OUTPUT_DIR=${2:-"output/01_qc"}
THREADS=${3:-4}

echo "=========================================="
echo "ATAC-seq FastQC 质量控制分析"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 FastQC 是否安装
if ! command -v fastqc &> /dev/null; then
    echo "错误: FastQC 未安装，请先安装 FastQC"
    exit 1
fi

# 运行 FastQC
echo "开始 FastQC 分析..."
fastqc -t ${THREADS} -o ${OUTPUT_DIR} ${INPUT_DIR}/*.fastq.gz ${INPUT_DIR}/*.fq.gz 2>/dev/null || \
fastqc -t ${THREADS} -o ${OUTPUT_DIR} ${INPUT_DIR}/*.fastq ${INPUT_DIR}/*.fq 2>/dev/null

echo "FastQC 分析完成!"
echo "结果保存在: ${OUTPUT_DIR}"

# 运行 MultiQC 汇总（如果安装）
if command -v multiqc &> /dev/null; then
    echo "正在生成 MultiQC 报告..."
    multiqc ${OUTPUT_DIR} -o ${OUTPUT_DIR}/multiqc
    echo "MultiQC 报告保存在: ${OUTPUT_DIR}/multiqc"
fi

# ATAC-seq 特异性 QC 检查
echo ""
echo "ATAC-seq 特异性检查:"
echo "- 检查是否有转座酶接头序列"
echo "- 检查片段大小分布（需要比对后分析）"
echo "- 建议运行 deepTools 的 plotFingerprint 评估信噪比"

echo "=========================================="
echo "质量控制分析完成!"
echo "=========================================="
