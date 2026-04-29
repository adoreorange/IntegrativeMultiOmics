#!/bin/bash

# ==============================
# ATAC-seq Motif Analysis with HOMER
# ==============================

# 1. 安装HOMER（若未安装，取消注释并运行）
# git clone https://github.com/kundajelab/homer.git
# cd homer
# make install

# 2. 参数设置
PEAK_FILE="atac_peaks.bed"  # ATAC-seq peak文件（BED格式，需包含染色体、起始、结束位置）
GENOME="mm10"              # 基因组版本（如mm10、hg38，需与peak文件匹配）
OUTPUT_DIR="homer_motif_results"  # 输出目录
THREADS=8                  # 并行线程数

# 3. 运行HOMER motif分析
findMotifsGenome.pl \
  ${PEAK_FILE} \
  ${GENOME} \
  ${OUTPUT_DIR} \
  -size given \          # 使用peak原始大小（而非固定长度）
  -p ${THREADS} \        # 并行计算
  -cache                 # 缓存基因组索引，加速后续分析

echo "HOMER motif分析完成！结果位于：${OUTPUT_DIR}"
echo "关键输出文件："
echo "  - knownResults/：已知motif匹配结果（如TF结合位点）"
echo "  - deNovo/：新发现的motif（logo图、序列）"
echo "  - homerResults.html：可视化报告（含motif logo和统计）"
