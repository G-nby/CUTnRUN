# ==============================================================================
#                      模块1: 数据预处理与清洗
# ==============================================================================
message("\n--- [模块1] 正在执行：数据预处理与清洗 ---")

# --- 检查是否需要跳过此步骤 ---
if (file.exists(cache_files$clean_counts)) {
  message("--> [跳过] 缓存文件 `01_clean_counts.rds` 已存在。")
} else {
  message("--> [运行] 正在生成干净的计数矩阵...")
  
  # 创建目录
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 加载所有库 (为独立运行)
  suppressPackageStartupMessages({
    library(tidyverse); library(GenomicRanges); library(readr)
  })
  
  # ==============================================================================
  # ## 第 3 节：数据加载与清洗 ##
  # ==============================================================================
  message("=== 第 3 节：开始加载与清洗数据 ===")
  
  # --- 自动创建输出目录 ---
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures/Correlation_Family_Level"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "figures/Correlation_repName_Level"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "figures/Correlation_Gene_Level"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "figures/Boxplots"), showWarnings = FALSE, recursive = TRUE)
  
  # --- 加载并初步处理原始计数 ---
  counts <- read_tsv(counts_file, skip = 1, col_names = TRUE, show_col_types = FALSE)
  counts <- counts[!is.na(counts$Start) & !is.na(counts$End),]
  counts_df <- counts %>% select(Geneid, 7:last_col())
  
  # --- 利用基因组黑名单进行数据清洗 ---
  message("--- 正在使用基因组黑名单过滤特征... ---")
  saf_annotations <- read_tsv(saf_file_path, col_names = c("GeneID", "Chr", "Start", "End", "Strand"), show_col_types = FALSE)
  saf_gr <- makeGRangesFromDataFrame(saf_annotations, keep.extra.columns = TRUE)
  
  denylist_df <- read_tsv(deny_list_file, col_names = c('Chr', 'Start', 'End'), show_col_types = FALSE)
  denylist_gr <- makeGRangesFromDataFrame(denylist_df)
  
  overlaps <- overlapsAny(saf_gr, denylist_gr)
  clean_geneids <- saf_gr[!overlaps]$GeneID
  
  counts_cleaned <- counts_df %>% filter(Geneid %in% clean_geneids)
  message(paste("--- 清洗完成，剩余", nrow(counts_cleaned), "个有效特征 (基因/TE) ---"))
  
  # --- 释放内存 ---
  rm(counts, counts_df, saf_annotations, saf_gr, denylist_df, denylist_gr, overlaps, clean_geneids)
  invisible(gc())
  message("--- 原始数据加载与清洗完成 ---\n")
  saveRDS(counts_cleaned, file = cache_files$clean_counts)
  message("--> [完成] 已保存清洗后的计数矩阵到 `01_clean_counts.rds`。")
}