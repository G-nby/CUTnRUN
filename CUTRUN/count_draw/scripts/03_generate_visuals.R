# ==============================================================================
#                      模块3: 生成所有可视化图表
# ==============================================================================
message("\n--- [模块3] 正在执行：生成所有可视化图表 ---")
message("--> [运行] 无论是否存在缓存，此模块总是重新生成图表。")

# 加载所有库
#suppressPackageStartupMessages({
#  library(tidyverse); library(ggplot2); library(ggrepel)
#  library(viridis); library(ggpubr); library(circlize)
#})

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel)
  library(viridis); library(circlize)
})

# --- 加载上一步的结果 ---
plot_data <- readRDS(cache_files$plot_data)

# --- 创建输出目录 ---
dir.create(file.path(figures_dir, "Correlation_Family_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_repName_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_Gene_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Boxplots"), recursive = TRUE, showWarnings = FALSE)
# ==============================================================================
# ## 第 6 节：定义核心绘图函数 ##
# ==============================================================================
message("=== 第 6 节：定义核心绘图函数 `generate_correlation_plot` ===")

#' @title 生成富集相关性散点图
#' @description 比较两个样本在基因或TE上的富集程度，并标注高差异目标
#' @param ... (参数与前一版相同)
#' @return 返回一个包含图中被标注数据点信息的数据框
generate_correlation_plot <- function(data, x_name, y_name, feature_type, te_classes = NULL,
                                      grouping_vars = NULL, label_var = NULL,
                                      n_labels, output_dir_path) {
  
  message(paste("    -> 正在绘制:", x_name, "vs", y_name, "| 类型:", feature_type))
  
  x_col <- sym(x_name)
  y_col <- sym(y_name)
  
  base_data <- data %>% filter(featureType == feature_type, sample %in% c(x_name, y_name))
  
  # --- TE 和 Gene 的数据处理逻辑分支 ---
  if (feature_type == "TE") {
    plot_df_base <- base_data %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarise(
        !!x_col := mean(lfc_over_igg[sample == x_name], na.rm = TRUE),
        !!y_col := mean(lfc_over_igg[sample == y_name], na.rm = TRUE),
        n_elements = n() / 2, mean_div = mean(milliDiv, na.rm = TRUE), .groups = "drop"
      ) %>%
      filter(Class %in% te_classes, n_elements > 10)
    
    if(nrow(plot_df_base) < n_labels) {
      message(paste("       警告: TE数据点不足 (", nrow(plot_df_base), ")，跳过绘图。"))
      return(NULL)
    }
    
    plot_df <- plot_df_base %>%
      mutate(age_group = cut(mean_div, breaks = quantile(mean_div, probs = seq(0, 1, 0.33), na.rm = TRUE),
                             labels = c("年轻", "中等", "古老"), include.lowest = TRUE))
    point_aes <- aes(color = Family, size = log10(n_elements), shape = age_group)
    size_lab <- "log10(元件数量)"; shape_lab <- "TE年龄组"; color_lab <- "TE家族"
    
  } else { # Gene
    plot_df <- base_data %>%
      select(sample, repName, Family, lfc_over_igg) %>%
      pivot_wider(names_from = sample, values_from = lfc_over_igg, values_fn = mean) %>%
      filter(!is.na(.data[[x_name]]) & !is.na(.data[[y_name]]))
    
    if(nrow(plot_df) < n_labels) {
      message(paste("       警告: Gene数据点不足 (", nrow(plot_df), ")，跳过绘图。"))
      return(NULL)
    }
    point_aes <- aes(color = Family)
    label_var <- "repName"
    size_lab <- NULL; shape_lab <- NULL; color_lab <- "基因生物类型"
  }
  
  # --- 筛选需要标注的数据点 (离原点最远的点) ---
  label_data <- plot_df %>%
    mutate(dist_from_origin = sqrt((!!x_col)^2 + (!!y_col)^2)) %>%
    slice_max(order_by = dist_from_origin, n = n_labels)
  
  cor_value <- cor(plot_df[[x_name]], plot_df[[y_name]], use = "complete.obs")
  p_value <- cor.test(plot_df[[x_name]], plot_df[[y_name]])$p.value
  label_text <- paste0("r = ", round(cor_value, 3), "\nP = ", signif(p_value, 3))

  # --- 核心绘图代码 ---
  p <- ggplot(data = plot_df, aes(x = !!x_col, y = !!y_col)) +
    geom_point(point_aes, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    #stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = label_text, size = 5) +
    geom_text_repel(
      data = label_data, aes(label = .data[[label_var]]),
      size = 3.5, box.padding = 0.5, min.segment.length = 0,
      segment.color = 'grey50', max.overlaps = 50, force = 3
    ) +
    scale_color_viridis_d(option = "turbo") +
    theme_bw(base_size = 14) +
    labs(
      title = paste(feature_type, "富集相关性:", x_name, "vs", y_name),
      subtitle = if(feature_type == "TE") paste("数据点按", label_var, "分组") else "每个数据点代表一个基因",
      x = paste(x_name, "的富集值 (log2FC vs IgG)"),
      y = paste(y_name, "的富集值 (log2FC vs IgG)"),
      size = size_lab, shape = shape_lab, color = color_lab
    )
  
  if (feature_type == "TE") { p <- p + scale_size_continuous(range = c(2, 10)) }
  
  # --- 保存图像 ---
  file_name <- paste0("Correlation_", feature_type, "_", x_name, "_vs_", y_name, ".png")
  ggsave(file.path(output_dir_path, file_name), p, width = 12, height = 10, dpi = 300, bg = "white")
  
  return(label_data)
}
message("--- 绘图函数定义完成 ---\n")


# ==============================================================================
# ## 第 7 节：执行批量相关性分析 ##
# ==============================================================================
message("=== 第 7 节：开始批量生成相关性分析图表 ===")

sample_list <- unique(plot_data$sample)
sample_pairs <- combn(sample_list, 2, simplify = FALSE)
all_significant_hits <- list()

for (pair in sample_pairs) {
  x_name <- pair[1]
  y_name <- pair[2]
  
  # --- a. TE 家族 (Family) 水平 ---
  family_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "TE",
    te_classes = TE_CLASSES_OF_INTEREST, grouping_vars = c("Family", "Class"),
    label_var = "Family", n_labels = N_LABELS_TE_FAMILY,
    output_dir_path = file.path(figures_dir, "Correlation_Family_Level")
  )
  if (!is.null(family_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      family_hits %>% mutate(feature_type = "TE", analysis_level = "Family", comparison = paste(x_name, "vs", y_name))
    ))
  }
  
  # --- b. TE 元件 (repName) 水平 ---
  repname_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "TE",
    te_classes = TE_CLASSES_OF_INTEREST, grouping_vars = c("repName", "Family", "Class"),
    label_var = "repName", n_labels = N_LABELS_TE_REPNAME,
    output_dir_path = file.path(figures_dir, "Correlation_repName_Level")
  )
  if (!is.null(repname_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      repname_hits %>% mutate(feature_type = "TE", analysis_level = "repName", comparison = paste(x_name, "vs", y_name))
    ))
  }
  
  # --- c. 基因 (Gene) 水平 ---
  gene_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "Gene",
    n_labels = N_LABELS_GENES,
    output_dir_path = file.path(figures_dir, "Correlation_Gene_Level")
  )
  if (!is.null(gene_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      gene_hits %>% mutate(feature_type = "Gene", analysis_level = "Gene", comparison = paste(x_name, "vs", y_name))
    ))
  }
}

# --- 汇总并保存所有图中被标注的"显著"目标 ---
if (length(all_significant_hits) > 0) {
  final_summary_df <- bind_rows(all_significant_hits)
  write_csv(final_summary_df, file.path(figures_dir, "significant_hits_summary.csv"))
  message("--- 所有相关性图表生成完毕，并已保存 `significant_hits_summary.csv` ---\n")
} else {
  message("--- 未能生成任何相关性图表或显著目标 ---\n")
}


# ==============================================================================
# ## 第 8 节：特定目标筛选与箱线图可视化 ##
# ==============================================================================
message("=== 第 8 节：开始根据指定条件筛选目标并绘制箱线图 ===")

# --- 8.1 TE筛选与绘图 ---
s1 <- SAMPLES_FOR_FILTERING[1]
s2 <- SAMPLES_FOR_FILTERING[2]
s3 <- SAMPLES_FOR_FILTERING[3]

te_df <- plot_data %>%
  filter(featureType == "TE", sample %in% SAMPLES_FOR_FILTERING) %>%
  group_by(repName, Family, Class) %>%
  summarise(
    !!sym(s1) := mean(lfc_over_igg[sample == s1], na.rm = TRUE),
    !!sym(s2) := mean(lfc_over_igg[sample == s2], na.rm = TRUE),
    !!sym(s3) := mean(lfc_over_igg[sample == s3], na.rm = TRUE),
    n_elements = n() / 3, .groups = "drop"
  ) %>%
  filter(n_elements > 10)

te_df_filtered <- te_df %>%
  filter(
    .data[[s1]] > median(te_df[[s1]], na.rm = TRUE),
    .data[[s3]] > median(te_df[[s3]], na.rm = TRUE),
    .data[[s2]] < median(te_df[[s2]], na.rm = TRUE)
  )

if (nrow(te_df_filtered) > 0) {
  filtered_te <- te_df_filtered %>%
    pivot_longer(cols = all_of(SAMPLES_FOR_FILTERING), names_to = "Sample", values_to = "lfc")
  
  p_family <- ggplot(filtered_te, aes(x = Sample, y = lfc, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(color = repName), width = 0.2, size = 2, alpha = 0.8) +
    facet_wrap(~Family, scales = "free_y") +
    labs(
      title = paste("在", s1, "和", s3, "中高富集，在", s2, "中低富集的TE"),
      x = "样本", y = "富集值 (log2FC vs IgG)", color = "TE 元件"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(file.path(figures_dir, "Boxplots/TE_Family_filtered_boxplot.png"), p_family, width = 16, height = 12, dpi = 300, bg = "white")
  message("--- 已生成筛选后的TE箱线图 (按家族划分) ---")
  
  # --- 图 8.1b: 按 TE 元件 (repName) 划分的点图 (新增加) ---
  p_repname <- ggplot(filtered_te, aes(x = Sample, y = lfc, color = Sample)) +
    # 使用点图展示每个样本的精确LFC值
    geom_point(size = 4, alpha = 0.8) +
    # 用线连接同一个repName在不同样本中的点，以展示变化趋势
    geom_line(aes(group = repName), color = "grey60", linetype = "dashed", alpha = 0.7) +
    # 每个面板(facet)代表一个独立的repName
    facet_wrap(~repName, scales = "free_y") +
    labs(
      title = paste("在", s1, "和", s3, "中高富集，在", s2, "中低富集的TE"),
      subtitle = "按TE元件(repName)独立展示",
      x = "样本", y = "富集值 (log2FC vs IgG)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none", # X轴和颜色信息重复，移除图例
      strip.background = element_rect(fill = "lightblue", color = "black"), # 美化分面标题
      strip.text = element_text(face = "bold", size = 10)
    )
  
  ggsave(
    file.path(figures_dir, "Boxplots/TE_repName_filtered_dotplot.png"),
    p_repname,
    width = 18,
    height = 18,
    dpi = 300,
    bg = "white"
  )
  message("--- 已生成筛选后的TE点图 (按元件repName划分) ---")
  
} else {
  message("--- 未找到符合筛选条件的TE ---")
}

# --- 8.2 基因筛选与绘图 ---
gene_df <- plot_data %>%
  filter(featureType == "Gene") %>%
  select(sample, repName, lfc_over_igg) %>%
  pivot_wider(names_from = sample, values_from = lfc_over_igg, values_fn = mean) %>%
  filter(if_all(all_of(SAMPLES_FOR_FILTERING), ~ !is.na(.)))

gene_df_filtered <- gene_df %>%
  filter(
    .data[[s1]] > median(gene_df[[s1]], na.rm = TRUE),
    .data[[s3]] > median(gene_df[[s3]], na.rm = TRUE),
    .data[[s2]] < median(gene_df[[s2]], na.rm = TRUE)
  ) %>%
  filter(.data[[s1]] > 0.7, .data[[s3]] > 0.7, .data[[s2]] < -0.2) # 更严格的阈值

if (nrow(gene_df_filtered) > 0) {
  samples_to_plot <- intersect(SAMPLES_FOR_GENE_BOXPLOT, colnames(gene_df_filtered))
  filtered_gene <- gene_df_filtered %>%
    pivot_longer(cols = all_of(samples_to_plot), names_to = "Sample", values_to = "lfc") %>%
    mutate(Sample = factor(Sample, levels = SAMPLES_FOR_GENE_BOXPLOT)) # 保持样本顺序
  
  p_gene <- ggplot(filtered_gene, aes(x = Sample, y = lfc, color = Sample)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_line(aes(group = repName), color = "grey80", alpha = 0.7) +
    facet_wrap(~repName, scales = "free_y") +
    labs(
      title = paste("在", s1, "和", s3, "中高富集，在", s2, "中低富集的基因"),
      x = "样本", y = "富集值 (log2FC vs IgG)"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(figures_dir, "Boxplots/Gene_filtered_dotplot.png"), p_gene, width = 14, height = 14, dpi = 300, bg = "white")
  message("--- 已生成筛选后的基因点图 ---")
} else {
  message("--- 未找到符合筛选条件的基因 ---")
}

message("--- 第 8 节分析完成 ---\n")

# ==============================================================================
# ## 第 9 节：分析结束 ##
# ==============================================================================
message("==========================================================")
message("          分析流程全部顺利完成！")
message(paste("          所有结果已保存至:", file.path(getwd(), output_dir)))
message("==========================================================")
