library(MOFA2)

# 读取数据
setwd('/home/hyf/rrbs')
gene <- read.xlsx("./data/gene_list.xlsx")
model <- run_mofa(
  gene,
  factors = 10,
  convergence_mode = "slow"
)
