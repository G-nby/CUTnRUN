message("\n--- [模块2] 正在执行：标准化与富集计算 ---")

# --- 检查是否需要跳过此步骤 ---
if (file.exists(cache_files$plot_data)) {
  message("--> [跳过] 缓存文件 `02_plot_data.rds` 已存在。")
} else {
  message("--> [运行] 正在进行标准化并计算LFC...")
  
  # 加载所有库
  suppressPackageStartupMessages({
    library(tidyverse); library(DESeq2); library(stringr)
  })
  # --- 加载上一步的结果 ---
  counts_cleaned <- readRDS(cache_files$clean_counts)
  
  # ==============================================================================
  # ## 第 4 节：计数标准化 ##
  # ==============================================================================
  message("=== 第 4 节：开始使用 DESeq2 进行计数标准化 ===")
  
  # --- 准备矩阵与样本信息 ---
  counts_mat <- counts_cleaned %>%
    column_to_rownames(var = "Geneid") %>%
    as.matrix()
  rm(counts_cleaned)
  #sample_info <- tibble(sample_full = colnames(counts_mat)) %>%
  #  mutate(sample_id = str_remove(sample_full, "_clean\\.bam$")) # 清理样本名
  sample_info <- tibble(sample_full = colnames(counts_mat)) %>%
     mutate(sample_id = basename(sample_full) %>% str_remove("_clean\\.bam$"))
  
  colnames(counts_mat) <- sample_info$sample_id
  
  colData <- data.frame(row.names = sample_info$sample_id)
  stopifnot(all(rownames(colData) == colnames(counts_mat))) # 一致性检查
  
  # --- 使用DESeq2估算大小因子并进行标准化 ---
  # design = ~1 表示我们不进行差异分析，仅利用DESeq2优秀的标准化算法。
  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = colData, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  counts_norm <- counts(dds, normalized = TRUE)
  
  
  # --- 根据设定的阈值过滤低信号特征 ---
  counts_norm_filtered <- counts_norm[rowSums(counts_norm) > MIN_TOTAL_NORMALIZED_COUNTS, ]
  
  rm(dds, counts_mat, counts_cleaned, colData)
  invisible(gc())
  
  message(paste("--- 标准化完成，过滤后剩余", nrow(counts_norm_filtered), "个特征 ---\n"))
  
  
  # ==============================================================================
  # ## 第 5 节：计算富集值并整合注释 ##
  # ==============================================================================
  message("=== 第 5 节：计算 LFC (log2富集倍数) 并整合注释信息 ===")
  
  # --- 转换为长格式数据以便处理 ---
  long_counts_norm <- as.data.frame(counts_norm_filtered) %>%
    rownames_to_column(var = "geneID") %>%
    pivot_longer(cols = -geneID, names_to = "sample", values_to = "norm_count")
  
  # --- 计算相对于IgG的LFC ---
  igg_counts <- long_counts_norm %>%
    filter(str_ends(sample, "IgG")) %>%
    group_by(geneID) %>%
    summarise(igg_count = mean(norm_count, na.rm = TRUE)) # 如果有多个IgG，取平均值
  
  final_data <- long_counts_norm %>%
    filter(!str_ends(sample, "IgG")) %>%
    left_join(igg_counts, by = "geneID") %>%
    mutate(lfc_over_igg = log2((norm_count + 1) / (igg_count + 1))) # +1避免log2(0)
  
  rm(long_counts_norm, igg_counts, counts_norm_filtered, counts_norm)
  invisible(gc())
  
  # --- 加载并合并基因与TE的详细注释 ---
  message("--- 正在整合基因与TE的注释信息... ---")
  full_repeats <- read_table(repeat_annotations_file, skip = 1, col_names = c("GeneID", "repName", "Class", "Family", "milliDiv"), show_col_types = FALSE)
  full_genes <- read_table(gene_annotations_file, show_col_types = FALSE)
  
  gene_data_processed <- final_data %>%
    filter(str_starts(geneID, "ENSG")) %>%
    left_join(full_genes, by = c("geneID" = "GeneID")) %>%
    mutate(featureType = "Gene", Class = NA_character_, Family = GeneBiotype, repName = GeneName) %>%
    filter(!is.na(GeneName) & !is.na(GeneBiotype))
  
  te_data_processed <- final_data %>%
    filter(!str_starts(geneID, "ENSG")) %>%
    left_join(full_repeats, by = c("geneID" = "GeneID")) %>%
    mutate(featureType = "TE") %>%
    filter(!is.na(Family))
  
  # --- 构建最终用于绘图的数据框 ---
  plot_data <- bind_rows(gene_data_processed, te_data_processed) %>%
    select(geneID, sample, lfc_over_igg, featureType, Class, Family, repName, milliDiv)
  
  rm(final_data, full_repeats, full_genes, gene_data_processed, te_data_processed)
  invisible(gc())
  
  # --- 保存预处理好的绘图数据，方便后续检查或复用 ---
  write_csv(plot_data, file.path(output_dir, "processed_data_for_plotting.csv"))
  message("--- 数据整合完成，已生成 `processed_data_for_plotting.csv` ---\n")
  
  # --- 再次检查用户配置的样本名是否有效 ---
  all_samples_in_data <- unique(plot_data$sample)
  user_samples <- c(SAMPLES_FOR_FILTERING, SAMPLES_FOR_GENE_BOXPLOT)
  missing_samples <- setdiff(user_samples, all_samples_in_data)
  if (length(missing_samples) > 0) {
    stop(paste("错误：您在配置区指定的样本名不存在于数据中，请检查拼写:", paste(missing_samples, collapse=", ")))
  }
  
  saveRDS(plot_data, file = cache_files$plot_data)
  message("--> [完成] 已保存最终绘图数据到 `02_plot_data.rds`。")
}
