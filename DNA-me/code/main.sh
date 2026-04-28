#!/bin/bash
# DNA 甲基化分析流程主脚本
# 用法: bash main.sh [config_file] [start_step] [end_step]

set -e

echo "============================================================"
echo "       DNA Methylation 分析流程"
echo "============================================================"
echo ""

# 默认配置
CONFIG_FILE=${1:-"config/parameters.conf"}
START_STEP=${2:-1}
END_STEP=${3:-5}

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
GENOME_FOLDER=${GENOME_FOLDER:-"${REFERENCE_DIR}/Bismark_genome"}

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# 显示流程信息
show_info() {
    echo ""
    echo "============================================================"
    echo "分析步骤:"
    echo "  1. 甲基化提取 (Bismark methylation extractor)"
    echo "  2. 甲基化计算 (methylKit)"
    echo "  3. DMR 分析 (DSS)"
    echo "  4. 富集分析 (GO/KEGG)"
    echo "  5. 结果可视化"
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

    # Bismark 相关工具
    tools=("bismark_methylation_extractor" "bismark2report")

    for tool in "${tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            log_warn "$tool 未安装"
        else
            log_info "$tool 已安装"
        fi
    done

    log_info "检查 R 环境..."
    Rscript -e "suppressPackageStartupMessages({library(methylKit); library(DSS)})" 2>/dev/null && \
        log_info "R 包已正确安装" || \
        log_warn "部分 R 包未安装"

    echo ""
}

# 步骤 0: Bismark 比对 (如果需要)
step0_alignment() {
    log_info "=========================================="
    log_info "步骤 0: Bismark 比对 (可选)"
    log_info "=========================================="

    if [ ! -d "${OUTPUT_DIR}/00_alignment" ]; then
        log_info "运行 Bismark 比对..."

        mkdir -p ${OUTPUT_DIR}/00_alignment

        # 检查 Bismark 基因组索引
        if [ ! -d "${GENOME_FOLDER}" ]; then
            log_info "构建 Bismark 基因组索引..."
            bismark_genome_preparation --path_to_bowtie2 $(which bowtie2 | xargs dirname) \
                ${REFERENCE_DIR}
        fi

        # 运行 Bismark 比对
        for R1 in ${DATA_DIR}/*_R1*.fastq.gz; do
            if [ -f "$R1" ]; then
                SAMPLE=$(basename $R1 | sed 's/_R1.*//')
                R2=$(echo $R1 | sed 's/_R1/_R2/')

                if [ -f "$R2" ]; then
                    log_info "比对样本: ${SAMPLE}"
                    bismark --genome ${GENOME_FOLDER} \
                        -1 ${R1} -2 ${R2} \
                        -o ${OUTPUT_DIR}/00_alignment/${SAMPLE} \
                        --multicore ${THREADS}
                fi
            fi
        done
    else
        log_info "BAM 文件已存在，跳过比对步骤"
    fi
}

# 步骤 1: 甲基化提取
step1_preprocess() {
    log_info "=========================================="
    log_info "步骤 1: 甲基化提取 (Bismark)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/01_preprocess
    bash code/01_preprocess/bismark_methylation_extractor.sh \
        ${OUTPUT_DIR}/00_alignment \
        ${OUTPUT_DIR}/01_preprocess \
        ${THREADS}

    log_info "步骤 1 完成!"
    echo ""
}

# 步骤 2: 甲基化计算
step2_calculation() {
    log_info "=========================================="
    log_info "步骤 2: 甲基化计算 (methylKit)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/02_calculation
    Rscript code/02_calculation/calculate_methylation.R \
        ${OUTPUT_DIR}/01_preprocess \
        ${OUTPUT_DIR}/02_calculation \
        data/sample_info.txt

    log_info "步骤 2 完成!"
    echo ""
}

# 步骤 3: DMR 分析
step3_dmr() {
    log_info "=========================================="
    log_info "步骤 3: DMR 分析 (DSS)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/03_dmr_analysis
    Rscript code/03_dmr_analysis/dss_dmr_analysis.R \
        ${OUTPUT_DIR}/02_calculation/methylation_objects.RData \
        data/sample_info.txt \
        ${OUTPUT_DIR}/03_dmr_analysis

    log_info "步骤 3 完成!"
    echo ""
}

# 步骤 4: 富集分析
step4_enrichment() {
    log_info "=========================================="
    log_info "步骤 4: 富集分析 (GO/KEGG)"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/04_enrichment
    Rscript code/04_enrichment/enrichment_analysis.R \
        ${OUTPUT_DIR}/03_dmr_analysis/DMR_significant.csv \
        ${OUTPUT_DIR}/04_enrichment

    log_info "步骤 4 完成!"
    echo ""
}

# 步骤 5: 可视化
step5_visualization() {
    log_info "=========================================="
    log_info "步骤 5: 结果可视化"
    log_info "=========================================="

    mkdir -p ${OUTPUT_DIR}/05_visualization
    Rscript code/05_visualization/visualization.R \
        ${OUTPUT_DIR}/03_dmr_analysis/dmr_analysis_results.RData \
        ${OUTPUT_DIR}/05_visualization

    log_info "步骤 5 完成!"
    echo ""
}

# 运行分析流程
run_pipeline() {
    local start_time=$(date +%s)

    show_info
    check_dependencies

    # 检查是否需要比对
    if [ ! -d "${OUTPUT_DIR}/00_alignment" ]; then
        log_warn "未找到比对结果，请先运行 Bismark 比对"
        log_warn "或者手动运行: bash code/00_alignment/bismark_align.sh"
    fi

    for step in $(seq $START_STEP $END_STEP); do
        case $step in
            1) step1_preprocess ;;
            2) step2_calculation ;;
            3) step3_dmr ;;
            4) step4_enrichment ;;
            5) step5_visualization ;;
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
    log_info "DNA Methylation 分析流程完成!"
    echo "============================================================"
    echo ""
    echo "总耗时: ${hours}h ${minutes}m ${seconds}s"
    echo "结果目录: ${OUTPUT_DIR}"
    echo ""
    echo "主要输出文件:"
    echo "  - 甲基化提取: ${OUTPUT_DIR}/01_preprocess/"
    echo "  - 甲基化矩阵: ${OUTPUT_DIR}/02_calculation/"
    echo "  - DMR 结果: ${OUTPUT_DIR}/03_dmr_analysis/"
    echo "  - 富集分析: ${OUTPUT_DIR}/04_enrichment/"
    echo "  - 可视化图表: ${OUTPUT_DIR}/05_visualization/"
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
        echo "  bash main.sh config/my_config.conf 2 5  # 从步骤2运行到步骤5"
        echo ""
        echo "可用步骤:"
        echo "  1 - 甲基化提取 (Bismark)"
        echo "  2 - 甲基化计算 (methylKit)"
        echo "  3 - DMR 分析 (DSS)"
        echo "  4 - 富集分析 (GO/KEGG)"
        echo "  5 - 结果可视化"
        echo ""
        echo "注意: 需要先运行 Bismark 比对生成 BAM 文件"
        exit 1
    fi

    run_pipeline
}

main "$@"
