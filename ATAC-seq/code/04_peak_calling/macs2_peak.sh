#!/bin/bash
# MACS2 Peak Calling (ATAC-seq)
# 用法: bash macs2_peak.sh <bam_dir> <output_dir>

set -e

# 参数设置
BAM_DIR=${1:-"output/03_alignment"}
OUTPUT_DIR=${2:-"output/04_peak_calling"}
GENOME_SIZE=${3:-"hs"}  # hs: 人类, mm: 小鼠, dm: 果蝇
QVALUE=${4:-0.01}

echo "=========================================="
echo "ATAC-seq MACS2 Peak Calling"
echo "=========================================="
echo "BAM 目录: $BAM_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "基因组大小: $GENOME_SIZE"
echo "Q-value: $QVALUE"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 检查 MACS2 是否安装
if ! command -v macs2 &> /dev/null; then
    echo "错误: MACS2 未安装"
    echo "安装方法: pip install macs2"
    exit 1
fi

# 检查 bedtools 是否安装
if ! command -v bedtools &> /dev/null; then
    echo "错误: bedtools 未安装"
    exit 1
fi

# ATAC-seq 峰调用
echo "开始 Peak Calling..."
echo ""

for BAM in ${BAM_DIR}/*/*.final.bam; do
    if [ -f "$BAM" ]; then
        SAMPLE=$(basename $BAM .final.bam)
        echo "处理样本: ${SAMPLE}"

        # MACS2 callpeak for ATAC-seq
        # --nomodel: 不建模，ATAC-seq 无移位
        # --shift -100 --extsize 200: ATAC-seq 推荐参数
        # --keep-dup all: 保留所有 reads (已去重)

        macs2 callpeak \
            -t ${BAM} \
            -f BAMPE \
            -g ${GENOME_SIZE} \
            -n ${SAMPLE} \
            -q ${QVALUE} \
            --nomodel \
            --shift -100 \
            --extsize 200 \
            --keep-dup all \
            --call-summits \
            -B \
            --SPMR \
            -outdir ${OUTPUT_DIR}/${SAMPLE}

        echo "  完成: ${SAMPLE}"
        echo ""
    fi
done

# 合并所有样本的 peaks
echo "合并所有样本的 peaks..."
cat ${OUTPUT_DIR}/*/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUTPUT_DIR}/merged_peaks.bed

# 生成 peak 统计
echo ""
echo "Peak 统计:"
echo "样本名,Peak数量,平均峰长" > ${OUTPUT_DIR}/peak_summary.csv

for DIR in ${OUTPUT_DIR}/*/; do
    SAMPLE=$(basename $DIR)
    if [ -f "${DIR}${SAMPLE}_peaks.narrowPeak" ]; then
        PEAK_COUNT=$(wc -l < ${DIR}${SAMPLE}_peaks.narrowPeak)
        AVG_LEN=$(awk '{sum+=$3-$2; count++} END {print int(sum/count)}' ${DIR}${SAMPLE}_peaks.narrowPeak)
        echo "${SAMPLE},${PEAK_COUNT},${AVG_LEN}"
        echo "${SAMPLE},${PEAK_COUNT},${AVG_LEN}" >> ${OUTPUT_DIR}/peak_summary.csv
    fi
done

# 计算 FRiP 分数
echo ""
echo "计算 FRiP 分数..."
echo "样本名,Total_reads,Peak_reads,FRiP" > ${OUTPUT_DIR}/frip_scores.csv

for BAM in ${BAM_DIR}/*/*.final.bam; do
    if [ -f "$BAM" ]; then
        SAMPLE=$(basename $BAM .final.bam)
        TOTAL=$(samtools view -c ${BAM})

        if [ -f "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_peaks.narrowPeak" ]; then
            # 计算 peak 区域的 reads 数量
            bedtools intersect -a ${BAM} -b ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_peaks.narrowPeak -u | \
                samtools view -c > /tmp/peak_reads.txt
            PEAK_READS=$(cat /tmp/peak_reads.txt)
            FRIP=$(echo "scale=4; ${PEAK_READS}/${TOTAL}" | bc)
            echo "${SAMPLE},${TOTAL},${PEAK_READS},${FRIP}"
            echo "${SAMPLE},${TOTAL},${PEAK_READS},${FRIP}" >> ${OUTPUT_DIR}/frip_scores.csv
        fi
    fi
done

# 生成 Tn5 插入位点 bed 文件（用于后续分析）
echo ""
echo "生成 Tn5 插入位点文件..."
for BAM in ${BAM_DIR}/*/*.final.bam; do
    if [ -f "$BAM" ]; then
        SAMPLE=$(basename $BAM .final.bam)
        echo "  处理: ${SAMPLE}"

        # 提取插入位点
        samtools view -b -F 0x4 ${BAM} | \
            bedtools bamtobed -i - | \
            awk 'BEGIN{OFS="\t"} $6=="+" {print $1, $2, $2+1, $4, ".", "+"} $6=="-" {print $1, $3-1, $3, $4, ".", "-"}' \
            > ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_insertion_sites.bed
    fi
done

echo ""
echo "=========================================="
echo "MACS2 Peak Calling 完成!"
echo "结果保存在: ${OUTPUT_DIR}"
echo ""
echo "主要输出文件:"
echo "  - *_peaks.narrowPeak: 窄峰"
echo "  - *_peaks.broadPeak: 宽峰 (如果有)"
echo "  - merged_peaks.bed: 合并的峰"
echo "  - peak_summary.csv: 峰统计"
echo "  - frip_scores.csv: FRiP 分数"
echo "=========================================="
