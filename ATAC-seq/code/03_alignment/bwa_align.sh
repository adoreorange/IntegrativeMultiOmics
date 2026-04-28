#!/bin/bash
# BWA 比对 (ATAC-seq)
# 用法: bash bwa_align.sh <input_dir> <output_dir> <genome_index>

set -e

# 参数设置
INPUT_DIR=${1:-"output/02_preprocess"}
OUTPUT_DIR=${2:-"output/03_alignment"}
GENOME_INDEX=${3:-"reference/BWA_index/genome"}
THREADS=${4:-8}

echo "=========================================="
echo "ATAC-seq BWA 序列比对"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "基因组索引: $GENOME_INDEX"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 BWA 是否安装
if ! command -v bwa &> /dev/null; then
    echo "错误: BWA 未安装"
    exit 1
fi

# 检查 samtools 是否安装
if ! command -v samtools &> /dev/null; then
    echo "错误: samtools 未安装"
    exit 1
fi

# 检查基因组索引是否存在
if [ ! -f "${GENOME_INDEX}.bwt" ]; then
    echo "错误: 基因组索引不存在: ${GENOME_INDEX}"
    echo "请先运行基因组索引构建:"
    echo "  bwa index -p reference/BWA_index/genome reference/genome.fa"
    exit 1
fi

# 处理配对末端数据
for R1 in ${INPUT_DIR}/*_R1_paired.fastq.gz; do
    if [ -f "$R1" ]; then
        SAMPLE=$(basename $R1 | sed 's/_R1_paired.*//')
        R2="${INPUT_DIR}/${SAMPLE}_R2_paired.fastq.gz"

        if [ -f "$R2" ]; then
            echo "比对样本: ${SAMPLE}"

            mkdir -p ${OUTPUT_DIR}/${SAMPLE}

            # BWA mem 比对
            echo "  运行 BWA mem..."
            bwa mem -t ${THREADS} -M ${GENOME_INDEX} ${R1} ${R2} 2>/dev/null | \
                samtools view -@ ${THREADS} -Sb - > ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.bam

            # 排序
            echo "  排序 BAM 文件..."
            samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.bam
            rm ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.bam

            # 建立索引
            echo "  建立索引..."
            samtools index ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam

            # 标记重复序列
            echo "  标记重复序列..."
            samtools markdup -@ ${THREADS} \
                ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam \
                ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.dedup.bam

            samtools index ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.dedup.bam

            # ATAC-seq 特异性处理：移除线粒体 reads
            echo "  移除线粒体 reads..."
            samtools idxstats ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.dedup.bam | \
                cut -f1 | grep -v -E "^chrM|^MT|^M" | xargs samtools view -@ ${THREADS} -b \
                ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.dedup.bam > \
                ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.final.bam

            samtools index ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.final.bam

            # 统计比对结果
            echo "  统计比对结果..."
            TOTAL=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam)
            DEDUP=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.dedup.bam)
            FINAL=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.final.bam)
            MT_READS=$((DEDUP - FINAL))

            echo "    Total reads: ${TOTAL}"
            echo "    After dedup: ${DEDUP}"
            echo "    Mitochondrial reads removed: ${MT_READS}"
            echo "    Final reads: ${FINAL}"
        fi
    fi
done

# 计算 FRiP 分数 (Fraction of Reads in Peaks) - 需要在 peak calling 后更新
echo ""
echo "注: FRiP 分数将在 peak calling 后计算"

# 生成比对统计报告
echo ""
echo "生成比对统计报告..."
echo "Sample,Total,Dedup,MT_removed,Final" > ${OUTPUT_DIR}/alignment_summary.csv

for DIR in ${OUTPUT_DIR}/*/; do
    SAMPLE=$(basename $DIR)
    if [ -f "${DIR}${SAMPLE}.sorted.bam" ]; then
        TOTAL=$(samtools view -c ${DIR}${SAMPLE}.sorted.bam)
        DEDUP=$(samtools view -c ${DIR}${SAMPLE}.dedup.bam)
        FINAL=$(samtools view -c ${DIR}${SAMPLE}.final.bam)
        MT_REMOVED=$((DEDUP - FINAL))
        echo "${SAMPLE},${TOTAL},${DEDUP},${MT_REMOVED},${FINAL}" >> ${OUTPUT_DIR}/alignment_summary.csv
    fi
done

echo ""
echo "=========================================="
echo "BWA 比对完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo "=========================================="
