#!/bin/bash
# Trimmomatic 质量修剪和接头去除
# 用法: bash trimmomatic.sh <input_dir> <output_dir>

set -e

# 参数设置
INPUT_DIR=${1:-"data/raw"}
OUTPUT_DIR=${2:-"output/02_preprocess"}
THREADS=${3:-4}
ADAPTER_FILE=${4:-"adapters/TruSeq3-PE-2.fa"}

echo "=========================================="
echo "Trimmomatic 质量修剪"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 Trimmomatic 是否安装
if ! command -v trimmomatic &> /dev/null; then
    echo "错误: Trimmomatic 未安装"
    exit 1
fi

# 处理配对末端数据
for R1 in ${INPUT_DIR}/*_R1*.fastq.gz ${INPUT_DIR}/*_1*.fastq.gz; do
    if [ -f "$R1" ]; then
        # 获取样本名
        SAMPLE=$(basename $R1 | sed 's/_R1.*//;s/_1.*//')

        # 查找对应的 R2 文件
        R2=$(echo $R1 | sed 's/_R1/_R2/;s/_1/_2/')

        if [ -f "$R2" ]; then
            echo "处理样本: ${SAMPLE}"

            trimmomatic PE -threads ${THREADS} -phred33 \
                ${R1} ${R2} \
                ${OUTPUT_DIR}/${SAMPLE}_R1_paired.fastq.gz ${OUTPUT_DIR}/${SAMPLE}_R1_unpaired.fastq.gz \
                ${OUTPUT_DIR}/${SAMPLE}_R2_paired.fastq.gz ${OUTPUT_DIR}/${SAMPLE}_R2_unpaired.fastq.gz \
                ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 \
                LEADING:3 \
                TRAILING:3 \
                SLIDINGWINDOW:4:15 \
                MINLEN:36
        else
            echo "警告: 未找到 ${SAMPLE} 的 R2 文件，跳过"
        fi
    fi
done

# 处理单端数据
for FQ in ${INPUT_DIR}/*.fastq.gz ${INPUT_DIR}/*.fq.gz; do
    if [ -f "$FQ" ] && [[ "$FQ" != *"_R1"* ]] && [[ "$FQ" != *"_R2"* ]] && [[ "$FQ" != *"_1"* ]] && [[ "$FQ" != *"_2"* ]]; then
        SAMPLE=$(basename $FQ .fastq.gz | sed 's/.fq.gz$//')
        echo "处理单端样本: ${SAMPLE}"

        trimmomatic SE -threads ${THREADS} -phred33 \
            ${FQ} \
            ${OUTPUT_DIR}/${SAMPLE}_trimmed.fastq.gz \
            ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36
    fi
done

echo "=========================================="
echo "Trimmomatic 处理完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo "=========================================="
