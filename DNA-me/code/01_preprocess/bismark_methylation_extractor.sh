#!/bin/bash
# Bismark 甲基化提取
# 用法: bash bismark_methylation_extractor.sh <bam_dir> <output_dir>

set -e

# 参数设置
BAM_DIR=${1:-"output/00_alignment"}
OUTPUT_DIR=${2:-"output/01_preprocess"}
THREADS=${3:-8}
GENOME_FOLDER=${4:-"reference/Bismark_genome"}

echo "=========================================="
echo "DNA Methylation - Bismark 提取"
echo "=========================================="
echo "BAM 目录: $BAM_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo "参考基因组: $GENOME_FOLDER"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 bismark_methylation_extractor 是否安装
if ! command -v bismark_methylation_extractor &> /dev/null; then
    echo "错误: bismark_methylation_extractor 未安装"
    echo "安装方法: conda install -c bioconda bismark"
    exit 1
fi

# 收集所有 BAM 文件
BAM_FILES=$(find ${BAM_DIR} -name "*.bam" | grep -v ".bai")

if [ -z "$BAM_FILES" ]; then
    echo "错误: 未找到 BAM 文件"
    echo "请先运行 Bismark 比对:"
    echo "  bismark --genome ${GENOME_FOLDER} -1 reads_R1.fq.gz -2 reads_R2.fq.gz -o alignment_output"
    exit 1
fi

echo "找到以下 BAM 文件:"
echo "$BAM_FILES"
echo ""

# 运行 bismark_methylation_extractor
echo "开始提取甲基化信息..."

for BAM in $BAM_FILES; do
    SAMPLE=$(basename $BAM .bam | sed 's/_bismark_bt2//' | sed 's/_pe//' | sed 's/_se//')
    echo "处理样本: ${SAMPLE}"

    # 提取甲基化信息
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --CX \
        --output ${OUTPUT_DIR}/${SAMPLE} \
        --no_overlap \
        --comprehensive \
        --merge_non_CG \
        --multicore ${THREADS} \
        ${BAM}

    echo "完成: ${SAMPLE}"
    echo ""
done

# 生成汇总报告
echo "生成汇总报告..."

# 合并所有样本的 CpG 报告
echo "合并 CpG 报告..."
cat ${OUTPUT_DIR}/*/*_bismark_bt2_pe.deduplicated.CX_report.txt > ${OUTPUT_DIR}/all_samples_CX_report.txt 2>/dev/null || \
cat ${OUTPUT_DIR}/*/*_bismark_bt2.deduplicated.CX_report.txt > ${OUTPUT_DIR}/all_samples_CX_report.txt 2>/dev/null || \
echo "警告: 未找到 CX_report 文件"

# 统计甲基化率
echo ""
echo "甲基化统计:"
echo "Sample,Total_C,C_methylated,Methylation_rate" > ${OUTPUT_DIR}/methylation_summary.csv

for DIR in ${OUTPUT_DIR}/*/; do
    SAMPLE=$(basename $DIR)

    # 查找 mbias 文件
    MBIAS=$(find $DIR -name "*M-bias.txt" | head -1)

    if [ -f "$MBIAS" ]; then
        # 从 M-bias 文件提取统计信息
        TOTAL_C=$(grep -A 10 "CpG context" $MBIAS | grep -v "^$" | tail -n +2 | awk '{sum+=$3+$4} END {print sum}')
        METHYL_C=$(grep -A 10 "CpG context" $MBIAS | grep -v "^$" | tail -n +2 | awk '{sum+=$3} END {print sum}')

        if [ ! -z "$TOTAL_C" ] && [ "$TOTAL_C" -gt 0 ]; then
            RATE=$(echo "scale=4; ${METHYL_C}/${TOTAL_C}" | bc)
            echo "${SAMPLE},${TOTAL_C},${METHYL_C},${RATE}"
            echo "${SAMPLE},${TOTAL_C},${METHYL_C},${RATE}" >> ${OUTPUT_DIR}/methylation_summary.csv
        fi
    fi
done

# 生成 Bismark 汇总报告
if command -v bismark2report &> /dev/null; then
    echo ""
    echo "生成 Bismark2 报告..."
    bismark2report --dir ${OUTPUT_DIR}/reports
fi

# 生成 BED 文件用于可视化
echo ""
echo "生成 BED 文件..."

for FILE in ${OUTPUT_DIR}/*/*.bedGraph.gz; do
    if [ -f "$FILE" ]; then
        SAMPLE=$(basename $FILE | sed 's/.bedGraph.gz//')
        echo "  处理: ${SAMPLE}"

        # 解压并转换为 BED 格式
        zcat $FILE | awk 'BEGIN{OFS="\t"} {if($4 >= 0) print $1, $2, $3, $4, ".", "+"}' \
            > ${OUTPUT_DIR}/${SAMPLE}_methylation.bed
    fi
done

echo ""
echo "=========================================="
echo "Bismark 甲基化提取完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo ""
echo "主要输出文件:"
echo "  - CX_report.txt: 甲基化报告"
echo "  - bedGraph.gz: 甲基化位点文件"
echo "  - methylation_summary.csv: 甲基化统计"
echo "=========================================="
