#!/usr/bin/env Rscript

# usage：
# Rscript convert_rds_to_expr.R /path/to/RE_all_3_TPM.RDS /path/to/output/EXPR.csv

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_csv <- args[2]

library(readr)
library(tibble)

# 加载 RDS 文件
rds <- readRDS(input_rds)

# 提取 counts（TPM 表达矩阵）
expr <- rds$counts  # rows: TE names, columns: samples

# 加入 TE 名作为第一列
expr_out <- rownames_to_column(as.data.frame(expr), var = "TE")

# 写为 CSV 文件，逗号分隔，header 带有 TE 和样本名
write_csv(expr_out, output_csv)
