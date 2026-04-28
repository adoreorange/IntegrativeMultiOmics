#!/bin/bash
# featureCounts 定量
# 用法: bash featurecounts.sh <bam_dir> <output_dir> <gtf_file>

set -e

# 参数设置
BAM_DIR=${1:-"output/03_alignment"}
OUTPUT_DIR=${2:-"output/04_quantification"}
GTF_FILE=${3:-"reference/genes.gtf"}
THREADS=${4:-4}

echo "=========================================="
echo "featureCounts 基因定量"
echo "=========================================="
echo "BAM 目录: $BAM_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "GTF 文件: $GTF_FILE"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 featureCounts 是否安装
if ! command -v featureCounts &> /dev/null; then
    echo "错误: featureCounts 未安装"
    exit 1
fi

# 检查 GTF 文件是否存在
if [ ! -f "$GTF_FILE" ]; then
    echo "错误: GTF 文件不存在: ${GTF_FILE}"
    exit 1
fi

# 收集所有 BAM 文件
BAM_FILES=$(find ${BAM_DIR} -name "*.sorted.bam" -o -name "*.bam" | grep -v ".bai" | sort)

if [ -z "$BAM_FILES" ]; then
    echo "错误: 未找到 BAM 文件"
    exit 1
fi

echo "找到以下 BAM 文件:"
echo "$BAM_FILES"

# 运行 featureCounts
echo "开始 featureCounts 定量..."

featureCounts -T ${THREADS} \
    -a ${GTF_FILE} \
    -o ${OUTPUT_DIR}/gene_counts.txt \
    -g gene_id \
    -t exon \
    -s 0 \
    -p \
    --countReadPairs \
    -B \
    -C \
    ${BAM_FILES}

# 生成基因计数矩阵（更简洁格式）
echo "生成基因计数矩阵..."
cut -f1,7- ${OUTPUT_DIR}/gene_counts.txt | tail -n +2 > ${OUTPUT_DIR}/count_matrix.txt

# 生成样本汇总统计
echo "生成样本统计..."
head -n 30 ${OUTPUT_DIR}/gene_counts.txt.summary > ${OUTPUT_DIR}/featurecounts_summary.txt

echo "=========================================="
echo "featureCounts 定量完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo "=========================================="
