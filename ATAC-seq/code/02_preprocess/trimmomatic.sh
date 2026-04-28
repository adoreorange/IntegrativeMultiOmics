#!/bin/bash
# Trimmomatic 质量修剪和接头去除 (ATAC-seq)
# ATAC-seq 使用 Tn5 转座酶，需要去除 Nextera 接头
# 用法: bash trimmomatic.sh <input_dir> <output_dir>

set -e

# 参数设置
INPUT_DIR=${1:-"data/raw"}
OUTPUT_DIR=${2:-"output/02_preprocess"}
THREADS=${3:-4}
# ATAC-seq 使用 Nextera 接头
ADAPTER_FILE=${4:-"adapters/NexteraPE-PE.fa"}

echo "=========================================="
echo "ATAC-seq Trimmomatic 质量修剪"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo "接头文件: $ADAPTER_FILE (Nextera)"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 Trimmomatic 是否安装
if ! command -v trimmomatic &> /dev/null; then
    echo "错误: Trimmomatic 未安装"
    exit 1
fi

# 创建 Nextera 接头文件（如果不存在）
if [ ! -f "$ADAPTER_FILE" ]; then
    echo "警告: Nextera 接头文件不存在，正在创建..."
    mkdir -p $(dirname $ADAPTER_FILE)
    cat > $ADAPTER_FILE << 'EOF'
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
CTGTCTCTTATACACATCT
>Trans1_rc
AGATGTGTATAAGAGACAG
>Trans2
CTGTCTCTTATACACATCT
>Trans2_rc
AGATGTGTATAAGAGACAG
EOF
    echo "已创建 Nextera 接头文件: $ADAPTER_FILE"
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

# 统计修剪结果
echo ""
echo "统计修剪结果:"
for file in ${OUTPUT_DIR}/*_paired.fastq.gz; do
    if [ -f "$file" ]; then
        count=$(zcat $file | echo $(( $(wc -l) / 4 )))
        echo "  $(basename $file): ${count} reads"
    fi
done

echo "=========================================="
echo "Trimmomatic 处理完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo "=========================================="
