#!/bin/bash
# RNA-seq 分析流程主脚本
# 用法: bash main.sh <config_file>

set -e

echo "============================================================"
echo "       RNA-seq 分析流程"
echo "============================================================"
echo ""

# 默认配置
CONFIG_FILE=${1:-"config/parameters.conf"}
START_STEP=${2:-1}
END_STEP=${3:-7}

# 加载配置文件
if [ -f "$CONFIG_FILE" ]; then
    echo "加载配置文件: $CONFIG_FILE"
    source $CONFIG_FILE
else
    echo "警告: 配置文件 $CONFIG_FILE 不存在，使用默认参数"
fi

# 设置默认参数
DATA_DIR=${DATA_DIR:-"data/raw"}
OUTPUT_DIR=${OUTPUT_DIR:-"output"}
REFERENCE_DIR=${REFERENCE_DIR:-"reference"}
THREADS=${THREADS:-8}
GTF_FILE=${GTF_FILE:-"${REFERENCE_DIR}/genes.gtf"}
STAR_INDEX=${STAR_INDEX:-"${REFERENCE_DIR}/STAR_index"}
ADAPTER_FILE=${ADAPTER_FILE:-"adapters/TruSeq3-PE-2.fa"}

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 显示流程信息
show_info() {
    echo ""
    echo "============================================================"
    echo "分析步骤:"
    echo "  1. 质量控制 (FastQC)"
    echo "  2. 数据预处理 (Trimmomatic)"
    echo "  3. 序列比对 (STAR)"
    echo "  4. 基因定量 (featureCounts)"
    echo "  5. 差异表达分析 (DESeq2)"
    echo "  6. 富集分析 (GO/KEGG)"
    echo "  7. 结果可视化"
    echo "============================================================"
    echo ""
    echo "当前配置:"
    echo "  数据目录: ${DATA_DIR}"
    echo "  输出目录: ${OUTPUT_DIR}"
    echo "  参考基因组: ${REFERENCE_DIR}"
    echo "  线程数: ${THREADS}"
    echo ""
    echo "执行步骤: ${START_STEP} - ${END_STEP}"
    echo "============================================================"
    echo ""
}

# 检查依赖
check_dependencies() {
    log_info "检查依赖软件..."

    # 检查 shell 脚本需要的工具
    tools=("fastqc" "trimmomatic" "STAR" "featureCounts" "samtools" "multiqc")

    for tool in "${tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            log_warn "$tool 未安装"
        else
            log_info "$tool 已安装"
        fi
    done

    # 检查 R 包
    log_info "检查 R 环境..."
    Rscript -e "suppressPackageStartupMessages({library(DESeq2); library(clusterProfiler)})" 2>/dev/null && \
        log_info "R 包已正确安装" || \
        log_warn "部分 R 包未安装，请运行: Rscript scripts/install_packages.R"

    echo ""
}

# 步骤 1: 质量控制
step1_qc() {
    log_info "=========================================="
    log_info "步骤 1: 质量控制 (FastQC)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/01_qc
    bash code/01_qc/fastqc.sh ${DATA_DIR} ${OUTPUT_DIR}/01_qc ${THREADS}

    log_info "步骤 1 完成!"
    echo ""
}

# 步骤 2: 数据预处理
step2_preprocess() {
    log_info "=========================================="
    log_info "步骤 2: 数据预处理 (Trimmomatic)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/02_preprocess
    bash code/02_preprocess/trimmomatic.sh ${DATA_DIR} ${OUTPUT_DIR}/02_preprocess ${THREADS} ${ADAPTER_FILE}

    log_info "步骤 2 完成!"
    echo ""
}

# 步骤 3: 序列比对
step3_alignment() {
    log_info "=========================================="
    log_info "步骤 3: 序列比对 (STAR)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/03_alignment
    bash code/03_alignment/star_align.sh ${OUTPUT_DIR}/02_preprocess ${OUTPUT_DIR}/03_alignment ${STAR_INDEX} ${THREADS}

    log_info "步骤 3 完成!"
    echo ""
}

# 步骤 4: 基因定量
step4_quantification() {
    log_info "=========================================="
    log_info "步骤 4: 基因定量 (featureCounts)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/04_quantification
    bash code/04_quantification/featurecounts.sh ${OUTPUT_DIR}/03_alignment ${OUTPUT_DIR}/04_quantification ${GTF_FILE} ${THREADS}

    log_info "步骤 4 完成!"
    echo ""
}

# 步骤 5: 差异表达分析
step5_de_analysis() {
    log_info "=========================================="
    log_info "步骤 5: 差异表达分析 (DESeq2)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/05_de_analysis
    Rscript code/05_de_analysis/deseq2_analysis.R \
        ${OUTPUT_DIR}/04_quantification/count_matrix.txt \
        data/sample_info.txt \
        ${OUTPUT_DIR}/05_de_analysis

    log_info "步骤 5 完成!"
    echo ""
}

# 步骤 6: 富集分析
step6_enrichment() {
    log_info "=========================================="
    log_info "步骤 6: 富集分析 (GO/KEGG)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/06_enrichment
    Rscript code/06_enrichment/enrichment_analysis.R \
        ${OUTPUT_DIR}/05_de_analysis/DEG_significant.csv \
        ${OUTPUT_DIR}/06_enrichment

    log_info "步骤 6 完成!"
    echo ""
}

# 步骤 7: 可视化
step7_visualization() {
    log_info "=========================================="
    log_info "步骤 7: 结果可视化"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/07_visualization
    Rscript code/07_visualization/visualization.R \
        ${OUTPUT_DIR}/05_de_analysis/deseq2_results.RData \
        ${OUTPUT_DIR}/07_visualization

    log_info "步骤 7 完成!"
    echo ""
}

# 运行分析流程
run_pipeline() {
    local start_time=$(date +%s)

    show_info
    check_dependencies

    for step in $(seq $START_STEP $END_STEP); do
        case $step in
            1) step1_qc ;;
            2) step2_preprocess ;;
            3) step3_alignment ;;
            4) step4_quantification ;;
            5) step5_de_analysis ;;
            6) step6_enrichment ;;
            7) step7_visualization ;;
            *) log_error "无效的步骤: $step" ;;
        esac
    done

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))

    echo ""
    echo "============================================================"
    log_info "RNA-seq 分析流程完成!"
    echo "============================================================"
    echo ""
    echo "总耗时: ${hours}h ${minutes}m ${seconds}s"
    echo "结果目录: ${OUTPUT_DIR}"
    echo ""
    echo "主要输出文件:"
    echo "  - 质量控制: ${OUTPUT_DIR}/01_qc/"
    echo "  - 比对结果: ${OUTPUT_DIR}/03_alignment/"
    echo "  - 表达矩阵: ${OUTPUT_DIR}/04_quantification/"
    echo "  - 差异基因: ${OUTPUT_DIR}/05_de_analysis/DEG_significant.csv"
    echo "  - 富集分析: ${OUTPUT_DIR}/06_enrichment/"
    echo "  - 可视化图表: ${OUTPUT_DIR}/07_visualization/"
    echo ""
    echo "============================================================"
}

# 主程序
main() {
    if [ "$#" -eq 0 ]; then
        echo "用法: bash main.sh [config_file] [start_step] [end_step]"
        echo ""
        echo "示例:"
        echo "  bash main.sh                          # 运行所有步骤"
        echo "  bash main.sh config/my_config.conf    # 使用指定配置文件"
        echo "  bash main.sh config/my_config.conf 3 7  # 从步骤3运行到步骤7"
        echo ""
        echo "可用步骤:"
        echo "  1 - 质量控制 (FastQC)"
        echo "  2 - 数据预处理 (Trimmomatic)"
        echo "  3 - 序列比对 (STAR)"
        echo "  4 - 基因定量 (featureCounts)"
        echo "  5 - 差异表达分析 (DESeq2)"
        echo "  6 - 富集分析 (GO/KEGG)"
        echo "  7 - 结果可视化"
        exit 1
    fi

    run_pipeline
}

main "$@"
