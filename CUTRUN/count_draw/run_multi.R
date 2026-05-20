# ==============================================================================
#                      主运行脚本 (run.R)
# ==============================================================================
# 此脚本将按顺序执行所有分析步骤。
# 首次运行时，它会完成所有计算并生成缓存文件。
# 后续运行时，已完成的步骤会自动跳过，只执行必要的计算。

#！！！将analysis文件夹所在的位置设置为工作文件夹：
#setwd("/Volumes/MCC/Zhaolab/CUT&RUN_analysis/analysis")

#source("config.R") # 加载所有路径和参数
config_file <- commandArgs(trailingOnly = TRUE)
if(length(config_file)) source(config_file) else source("config.R")

message("=== 第 1 节：开始环境设置与依赖包检查 ===")
# --- 检查并安装CRAN包 ---
#cran_packages <- c("tidyverse", "ggrepel", "viridis", "ggpubr", "devtools")
cran_packages <- c("tidyverse", "ggrepel", "viridis", "devtools")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("正在从CRAN安装:", pkg, "..."))
    install.packages(pkg)
  }
}

# --- 检查并安装Bioconductor包 ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_packages <- c("GenomicRanges", "ComplexHeatmap", "circlize", "DESeq2")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("正在从Bioconductor安装:", pkg, "..."))
    BiocManager::install(pkg)
  }
}

# --- 静默加载所有库 ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(viridis)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(DESeq2)
})
message("--- 所有依赖包已成功加载 ---\n")

# --- 2. 依次执行各个分析模块 ---
message("=========================================")
message("### 开始执行分析流程 ###")
message("=========================================")

# 执行第1步：数据预处理
#source("scripts/01_prepare_data.R")

# 执行第2步：标准化与富集计算
#source("scripts/02_normalize_and_process.R")

# 执行第3步：可视化
#source("scripts/03_generate_visuals.R")
#this_dir <- dirname(sys.frame(1)$ofile) 
#this_dir <- dirname(/Share/home/gby/CUTRUN/
#suppressPackageStartupMessages(library(R.utils))
#this_dir <- dirname(thisFile())
#this_dir <- dirname(normalizePath(thisFile()))
message(counts_file)
cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", grep("--file=", cmd_args, value = TRUE))
if (length(script_path) == 1) {
    this_dir <- dirname(normalizePath(script_path))
} 

#source(file.path(this_dir, "scripts", "01_prepare_data.R"))
#source(file.path(this_dir, "scripts", "02_normalize_and_process.R"))
#source(file.path(this_dir, "scripts", "03_generate_visuals.R"))
message("\n=========================================")
message("### 开始执行多批次并行分析流程 ###")
message("=========================================")

base_intermediate <- intermediate_dir
merged_plot_data_list <- list()

# 循环处理每个组的 featureCounts 结果
for (gname in names(MULTI_COUNTS_FILES)) {
    message(paste("\n>>> 正在独立处理组别:", gname, "<<<"))
    
    # 动态覆盖 01 和 02 脚本需要的全局变量
    counts_file <<- MULTI_COUNTS_FILES[[gname]]
    
    # 为每组创建独立的缓存目录，防止互相覆盖
    group_intermediate <<- file.path(base_intermediate, gname)
    dir.create(group_intermediate, recursive = TRUE, showWarnings = FALSE)
    
    cache_files <<- list(
      clean_counts = file.path(group_intermediate, "01_clean_counts.rds"),
      plot_data    = file.path(group_intermediate, "02_plot_data.rds")
    )
    
    # 执行原生 01 和 02 脚本
    source(file.path(this_dir, "scripts", "01_prepare_data.R"))
    source(file.path(this_dir, "scripts", "02_normalize_and_process.R"))
    
    # 读取跑完的 02 结果，加上组名后缀
    group_plot_data <- readRDS(cache_files$plot_data) %>%
        mutate(sample = paste0(sample, "_", gname))
    
    merged_plot_data_list[[gname]] <- group_plot_data
}

message("\n>>> 所有独立组别处理完毕，正在合并数据... <<<")
plot_data_final <- bind_rows(merged_plot_data_list)

# 将合并后的数据存入总的 intermediate 目录
merged_rds_path <- file.path(base_intermediate, "merged_plot_data.rds")
saveRDS(plot_data_final, merged_rds_path)

# 将 cache_files 的指针修改为合并后的数据，欺骗 03 脚本直接读取
cache_files <<- list(plot_data = merged_rds_path)

message(">>> 进入绘图模块 <<<")
# 运行修改后的多批次绘图脚本
source(file.path(this_dir, "scripts", "03_generate_visuals_multi.R"))

message("=========================================")
message("### 多批次分析流程全部执行完毕！ ###")
message("=========================================")

message("=========================================")
message("### 分析流程全部执行完毕！ ###")
message("=========================================")
