#!/bin/bash
# STAR 比对
# 用法: bash star_align.sh <input_dir> <output_dir> <genome_dir>

set -e

# 参数设置
INPUT_DIR=${1:-"output/02_preprocess"}
OUTPUT_DIR=${2:-"output/03_alignment"}
GENOME_DIR=${3:-"reference/STAR_index"}
THREADS=${4:-8}

echo "=========================================="
echo "STAR 序列比对"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "基因组索引: $GENOME_DIR"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 STAR 是否安装
if ! command -v STAR &> /dev/null; then
    echo "错误: STAR 未安装"
    exit 1
fi

# 检查基因组索引是否存在
if [ ! -d "$GENOME_DIR" ]; then
    echo "错误: 基因组索引目录不存在: ${GENOME_DIR}"
    echo "请先运行基因组索引构建"
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

            STAR --runThreadN ${THREADS} \
                --genomeDir ${GENOME_DIR} \
                --readFilesIn ${R1} ${R2} \
                --readFilesCommand zcat \
                --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}/ \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode GeneCounts \
                --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} LB:lib1 PL:ILLUMINA \
                --outFilterMultimapNmax 20 \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000

            # 重命名输出文件
            mv ${OUTPUT_DIR}/${SAMPLE}/Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam
            mv ${OUTPUT_DIR}/${SAMPLE}/ReadsPerGene.out.tab ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_gene_counts.tab 2>/dev/null || true

            # 建立索引
            samtools index ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam
        fi
    fi
done

# 汇总比对统计
echo "汇总比对统计信息..."
echo "Sample,Uniquely_mapped,Multimapped,Too_many_mismatches,Too_short,Other" > ${OUTPUT_DIR}/alignment_summary.csv

for DIR in ${OUTPUT_DIR}/*/; do
    SAMPLE=$(basename $DIR)
    if [ -f "${DIR}LogFinal.out" ]; then
        UNIQUE=$(grep "Uniquely mapped reads number" ${DIR}LogFinal.out | awk '{print $NF}')
        MULTI=$(grep "Number of reads mapped to multiple loci" ${DIR}LogFinal.out | awk '{print $NF}')
        MISMATCH=$(grep "Number of reads mapped to too many loci" ${DIR}LogFinal.out | awk '{print $NF}')
        SHORT=$(grep "% of reads unmapped: too short" ${DIR}LogFinal.out | awk '{print $NF}')
        OTHER=$(grep "% of reads unmapped: other" ${DIR}LogFinal.out | awk '{print $NF}')
        echo "${SAMPLE},${UNIQUE},${MULTI},${MISMATCH},${SHORT},${OTHER}" >> ${OUTPUT_DIR}/alignment_summary.csv
    fi
done

echo "=========================================="
echo "STAR 比对完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo "=========================================="
