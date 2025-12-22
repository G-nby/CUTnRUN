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
  library(viridis); library(circlize); library(patchwork)
})

# --- 加载上一步的结果 ---
plot_data <- readRDS(cache_files$plot_data)

# --- 创建输出目录 ---
#dir.create(file.path(figures_dir, "Correlations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_Family_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_repName_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_Gene_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Boxplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Heatmaps"), recursive = TRUE, showWarnings = FALSE)

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
      ) 
      #%>%
      #filter(Class %in% te_classes, n_elements > 10)
    if (!identical(te_classes, "all")) {
      plot_df_base <- plot_df_base %>% filter(Class %in% te_classes)
    }
    
    plot_df_base <- plot_df_base %>% filter(n_elements > 10)
    
    if(nrow(plot_df_base) < n_labels) {
      message(paste("       警告: TE数据点不足 (", nrow(plot_df_base), ")，跳过绘图。"))
      return(NULL)
    }
    
    plot_df <- plot_df_base %>%
      mutate(age_group = cut(mean_div, breaks = quantile(mean_div, probs = seq(0, 1, 0.33), na.rm = TRUE),
                             labels = c("young", "mid", "old"), include.lowest = TRUE))
    point_aes <- aes(color = Family, size = log10(n_elements), shape = age_group)
    size_lab <- "log10(element count)"; shape_lab <- "TE age group"; color_lab <- "TE family"
    
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
    size_lab <- NULL; shape_lab <- NULL; color_lab <- "gene biotype"
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
      title = paste(feature_type, "enrichment correlation:", x_name, "vs", y_name),
      subtitle = if(feature_type == "TE") paste("data point are grouped by", label_var) else "every data point represent a gene",
      x = paste(x_name, " enrichment (log2FC vs IgG)"),
      y = paste(y_name, " enrichment (log2FC vs IgG)"),
      size = size_lab, shape = shape_lab, color = color_lab
    )
  
  if (feature_type == "TE") { p <- p + scale_size_continuous(range = c(2, 10)) }
  
  # --- 保存图像 ---
  file_name <- paste0("Correlation_", feature_type, "_", x_name, "_vs_", y_name, ".png")
  ggsave(file.path(output_dir_path, file_name), p, width = 12, height = 10, dpi = 300, bg = "white")
  
  return(label_data)
  #if (for_matrix) {
  #  p <- p +
  #    theme(
  #      legend.position = "none",
  #      plot.title = element_text(size = 10),
  #      axis.title = element_blank()
  #    )
  #}
  #return(list(data = label_data, plot = p))
  
}





generate_correlation_matrix_plot <- function(data, feature_type,
                                             te_classes=NULL, grouping_vars=NULL,
                                             label_var=NULL, n_labels=10,
                                             output_file, matrix_samples=NULL) {
  
  message("生成整体矩阵图...")
  
  if (is.null(matrix_samples)) matrix_samples <- unique(data$sample)
  sample_pairs <- combn(matrix_samples, 2, simplify = FALSE)
  
  all_plot_data <- list()
  all_label_data <- list()
  
  for (pair in sample_pairs) {
    x_name <- pair[1]; y_name <- pair[2]
    x_col <- sym(x_name); y_col <- sym(y_name)
    
    base_data <- data %>% filter(featureType == feature_type, sample %in% c(x_name, y_name))
    
    if (feature_type=="TE") {
      plot_df <- base_data %>%
        group_by(across(all_of(grouping_vars))) %>%
        summarise(
          x_val = mean(lfc_over_igg[sample==x_name], na.rm=TRUE),
          y_val = mean(lfc_over_igg[sample==y_name], na.rm=TRUE),
          n_elements = n()/2,
          mean_div = mean(milliDiv, na.rm=TRUE),
          .groups="drop"
        )
      if (!identical(te_classes,"all")) plot_df <- plot_df %>% filter(Class %in% te_classes)
      plot_df <- plot_df %>% filter(n_elements>10)
      
      if (nrow(plot_df)<1) next
      
      plot_df <- plot_df %>%
        mutate(age_group=cut(mean_div, breaks=quantile(mean_div, probs=seq(0,1,0.33), na.rm=TRUE),
                             labels=c("young","mid","old"), include.lowest=TRUE))
    } else { # Gene 
      plot_df <- base_data %>%
      select(sample, repName, Family, lfc_over_igg) %>%
      pivot_wider(names_from = sample, values_from = lfc_over_igg, values_fn = mean) %>%
      filter(!is.na(.data[[x_name]]) & !is.na(.data[[y_name]])) %>%
      rename_with(~c("x_val","y_val"), .cols = c(x_name, y_name))
      if (nrow(plot_df)<1) next
    }
    
    
    plot_df <- plot_df %>% mutate(pair=paste(x_name, y_name, sep="_"))
    
    label_data <- plot_df %>%
      mutate(dist_from_origin=sqrt(x_val^2+y_val^2)) %>%
      slice_max(order_by=dist_from_origin, n=n_labels) %>%
      mutate(pair=paste(x_name,y_name,sep="_"))
    
    all_plot_data[[length(all_plot_data)+1]] <- plot_df
    all_label_data[[length(all_label_data)+1]] <- label_data
  }
  
  combined_data <- bind_rows(all_plot_data)
  combined_labels <- bind_rows(all_label_data)
  
  # 如果没有有效数据，直接退出
  if (nrow(combined_data)==0) {
    message("无有效数据，跳过绘图: ", output_file)
    return(NULL)
  }
  
  # 为每行生成 row/col
  combined_data <- combined_data %>%
    rowwise() %>%
    mutate(row = factor(strsplit(pair,"_")[[1]][1], levels=matrix_samples),
           col = factor(strsplit(pair,"_")[[1]][2], levels=matrix_samples)) %>%
    ungroup()
  
  combined_labels <- combined_labels %>%
    rowwise() %>%
    mutate(row = factor(strsplit(pair,"_")[[1]][1], levels=matrix_samples),
           col = factor(strsplit(pair,"_")[[1]][2], levels=matrix_samples)) %>%
    ungroup()
  
  # --- 设置 aes ---
  if (feature_type=="TE") {
    geom_point_aes <- aes(color=Family, size=log10(n_elements), shape=age_group)
  } else {
    geom_point_aes <- aes(color=Family)
  }
  
  p <- ggplot(combined_data, aes(x=x_val, y=y_val)) +
    geom_point(geom_point_aes, alpha=0.7) +
    geom_vline(xintercept=0, linetype="dotted", color="grey50") +
    geom_hline(yintercept=0, linetype="dotted", color="grey50") +
    geom_text_repel(data=combined_labels,
                    aes(x=x_val, y=y_val, label=.data[[label_var]]), size=2.5) +
    facet_grid(row~col, scales="fixed") +
    scale_color_viridis_d(option="turbo") +
    theme_bw(base_size=14) +
    theme(strip.text=element_text(size=8),
          axis.title=element_blank(),
          legend.position="bottom")
  
  ggsave(output_file, p, width=20, height=20, bg="white")
  message("矩阵图保存完成: ", output_file)
}





message("--- 绘图函数定义完成 ---\n")


# ==============================================================================
# ## 第 7 节：执行批量相关性分析 ##
# ==============================================================================
message("=== 第 7 节：开始批量生成相关性分析图表 ===")

sample_list <- unique(plot_data$sample)
sample_pairs <- combn(sample_list, 2, simplify = FALSE)
all_significant_hits <- list()
all_plots_te_family <- list()

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


# 调用函数
generate_correlation_matrix_plot(
  data = plot_data,              
  feature_type = "TE",          
  te_classes = TE_CLASSES_OF_INTEREST, 
  grouping_vars = c("Family", "Class"), 
  label_var = "Family",        
  n_labels = N_LABELS_TE_FAMILY,
  output_file =  file.path(figures_dir, "Correlation_Family_Level/TE_Family_Correlation_Matrix.pdf"),    
  matrix_samples = unique(plot_data$sample) 
)

generate_correlation_matrix_plot(
  data = plot_data,
  feature_type = "TE",
  te_classes = TE_CLASSES_OF_INTEREST,
  grouping_vars = c("repName", "Family", "Class"),
  label_var = "repName",
  n_labels = N_LABELS_TE_REPNAME,
  output_file = file.path(figures_dir, "Correlation_repName_Level/TE_repName_Correlation_Matrix.pdf"),
  matrix_samples = unique(plot_data$sample)
)

generate_correlation_matrix_plot(
  data = plot_data,
  feature_type = "Gene",
  grouping_vars = c("repName"),
  label_var = "repName",
  n_labels = N_LABELS_GENES,
  output_file = file.path(figures_dir, "Correlation_Gene_Level/Gene_Correlation_Matrix.pdf"),
  matrix_samples = unique(plot_data$sample)
)



# ==============================================================================
# ## 第 8 节：特定目标筛选与箱线图可视化 ##
# ==============================================================================
#message("=== 第 8 节：绘制箱线图 ===")
message("=== 第 8 节：TE & Gene Heatmap 和 Boxplot 可视化 ===")

if(!exists("TE_repname_OI")) TE_repname_OI <- NULL

# --- TE 数据 ---
te_data <- plot_data %>%
  filter(featureType=="TE", Class %in% TE_CLASSES_OF_INTEREST)
if(!is.null(TE_repname_OI) && length(TE_repname_OI)>0 && any(TE_repname_OI != "")){
  te_data_rep <- te_data %>% filter(repName %in% TE_repname_OI)
} else {
  te_data_rep <- te_data
}
message(TE_repname_OI)

# --- Gene 数据 ---
gene_data <- plot_data %>%
  filter(featureType=="Gene")

# --- 函数：计算相关性矩阵 ---
compute_cor_matrix <- function(df, value_col="lfc_over_igg"){
  samples <- unique(df$sample)
  cor_mat <- matrix(NA, nrow=length(samples), ncol=length(samples),
                    dimnames=list(samples, samples))
  sample_pairs <- combn(samples,2,simplify=FALSE)
  for(pair in sample_pairs){
    x <- pair[1]; y <- pair[2]
    x_vals <- df %>% filter(sample==x) %>% pull({{value_col}})
    y_vals <- df %>% filter(sample==y) %>% pull({{value_col}})
    cor_mat[x,y] <- cor(x_vals, y_vals, use="complete.obs")
    cor_mat[y,x] <- cor_mat[x,y]
  }
  diag(cor_mat) <- 1
  as.data.frame(as.table(cor_mat)) %>%
    setNames(c("Sample1","Sample2","Correlation"))
}

# ================= TE 图 =================
if(nrow(te_data)>0){

  samples <- unique(te_data$sample)

  # --- 总体 Heatmap ---
  cor_long <- compute_cor_matrix(te_data)
  p_heatmap <- ggplot(cor_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=3) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    theme_minimal() +
    labs(title="TE Correlation Heatmap (All TE Classes)")
  ggsave(file.path(figures_dir,"Heatmaps","TE_All_Heatmap.png"), p_heatmap, width=8, height=6)

  # --- 总体 Boxplot ---
  p_box <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    theme_bw(base_size=14) +
    labs(title="TE log2FC distribution across samples", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","TE_All_Boxplot.png"), p_box, width=8, height=6)

  # --- 按 Family (Class facet) Boxplot ---
  p_box_family <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~Family, scales="free_y") +
    theme_bw(base_size=12) +
    labs(title="TE log2FC by Family", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_Family.png"), p_box_family, width=12, height=8)
  
  p_box_class <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~Class, scales="free_y") +
    theme_bw(base_size=12) +
    labs(title="TE log2FC by Class", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_Class.png"), p_box_class, width=12, height=8)

  # --- 按 repName Boxplot ---
  p_box_repname <- ggplot(te_data_rep, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~repName, scales="free_y") +
    theme_bw(base_size=10) +
    labs(title="TE log2FC by repName", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_repName.png"), p_box_repname, width=16, height=12)

  # --- 按 Class facet Heatmap ---
  cor_class_list <- list()
  for(cl in unique(te_data$Class)){
    df_cl <- te_data %>% filter(Class==cl)
    cor_cl <- compute_cor_matrix(df_cl)
    cor_cl$Class <- cl
    cor_class_list[[length(cor_class_list)+1]] <- cor_cl
  }
  cor_class_long <- bind_rows(cor_class_list)
  p_heatmap_class <- ggplot(cor_class_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=2.5) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    facet_wrap(~Class) +
    theme_minimal() +
    labs(title="TE Correlation Heatmap by Class")
  ggsave(file.path(figures_dir,"Heatmaps","TE_Heatmap_by_Class.png"), p_heatmap_class, width=12, height=10)
  
  cor_family_list <- list()
  for(fa in unique(te_data$Family)){
    df_fa <- te_data %>% filter(Family == fa)
    cor_fa <- compute_cor_matrix(df_fa)
    cor_fa$Family <- fa
    cor_family_list[[length(cor_family_list)+1]] <- cor_fa
  }
  cor_family_long <- bind_rows(cor_family_list)
  p_heatmap_family <- ggplot(cor_family_long, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.5) +
    scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
    facet_wrap(~Family) +
    theme_minimal() +
    labs(title = "TE Correlation Heatmap by Family")
  ggsave(file.path(figures_dir, "Heatmaps", "TE_Heatmap_by_Family.png"), p_heatmap_family, width = 12, height = 10)

  # --- 按 repName facet Heatmap ---
  cor_rep_list <- list()
  for(rep in unique(te_data_rep$repName)){
    df_rep <- te_data %>% filter(repName==rep)
    cor_rep <- compute_cor_matrix(df_rep)
    cor_rep$repName <- rep
    cor_rep_list[[length(cor_rep_list)+1]] <- cor_rep
  }
  cor_rep_long <- bind_rows(cor_rep_list)
  p_heatmap_rep <- ggplot(cor_rep_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=2.5) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    facet_wrap(~repName) +
    theme_minimal() +
    labs(title="TE Correlation Heatmap by repName")
  ggsave(file.path(figures_dir,"Heatmaps","TE_Heatmap_by_repName.png"), p_heatmap_rep, width=16, height=12)

  message("--- TE Heatmap 和 Boxplot 已完成 ---")
}

# ================= Gene 图 =================
if(nrow(gene_data)>0){

  samples <- unique(gene_data$sample)

  # --- Gene 总体 Heatmap ---
  cor_long_gene <- compute_cor_matrix(gene_data)
  p_heatmap_gene <- ggplot(cor_long_gene, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=3) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    theme_minimal() +
    labs(title="Gene Correlation Heatmap")
  ggsave(file.path(figures_dir,"Heatmaps","Gene_All_Heatmap.png"), p_heatmap_gene, width=8, height=6)

  # --- Gene 总体 Boxplot ---
  p_box_gene <- ggplot(gene_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    theme_bw(base_size=14) +
    labs(title="Gene log2FC distribution across samples", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","Gene_All_Boxplot.png"), p_box_gene, width=8, height=6)

  message("--- Gene Heatmap 和 Boxplot 已完成 ---")
}


message("--- 第 8 节分析完成 ---\n")

# ==============================================================================
# ## 第 9 节：分析结束 ##
# ==============================================================================
message("==========================================================")
message("          分析流程全部顺利完成！")
message(paste("          所有结果已保存至:", file.path(getwd(), output_dir)))
message("==========================================================")