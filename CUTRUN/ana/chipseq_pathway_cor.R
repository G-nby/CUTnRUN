#!/opt/anaconda3/envs/gogsea/bin/Rscript

#required_packages <- c("dplyr", "ggplot2", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

#for (pkg in required_packages) {
#  if (!requireNamespace(pkg, quietly = TRUE)) {
#    message(paste("Installing", pkg, "..."))
#    if (pkg %in% c("dplyr", "ggplot2")) install.packages(pkg)
#    else BiocManager::install(pkg)
#  }
#  library(pkg, character.only = TRUE, quietly = TRUE)
#}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("用法: Rscript script.R <input_dir> <threshold> [output_dir]")
}

input_dir <- args[1]
threshold <- as.numeric(args[2])
if (is.na(threshold)) stop("请输入有效的数字作为 threshold")

output_dir <- ifelse(length(args) >= 3, args[3], file.path(input_dir, "pathway_analysis"))

if (!dir.exists(output_dir)) dir.create(output_dir)


counts_path <- file.path(input_dir, "01_clean_counts.rds")
lfc_path <- file.path(input_dir, "02_plot_data.rds")

if (!file.exists(counts_path) | !file.exists(lfc_path)) {
  stop("no required files! please check the params")
}

data1 <- readRDS(counts_path) # 原始 count
data <- readRDS(lfc_path)    # LFC 数据

#head(data)
#head(data1)

igg_col <- grep("IgG", colnames(data1), value = TRUE)
if (length(igg_col) == 0) stop("在 data1 中未找到 IgG 列，请检查列名")

#print(class(data1))
#print(colnames(data1)[1:3]) # 查看前三列的名字


#keep_genes <- data1[rowSums(data1[,-1]) > 10, 1] %>% as.character()
id_vector <- as.character(data1[[1]])


numeric_cols <- data1[, sapply(data1, is.numeric)]
row_counts <- rowSums(numeric_cols)

keep_genes <- id_vector[row_counts > 10]
keep_genes <- trimws(keep_genes) 
keep_genes <- keep_genes[grepl("^ENSG", keep_genes)]
  
#keep_genes <- keep_genes[grepl("^ENSG", keep_genes)]
#print(keep_genes)

message(paste("Total valid ENSG genes found in data1:", length(keep_genes)))
#print(head(keep_genes))
#head(keep_genes[0:5])
#str(keep_genes)
#type(keep_genes)

all_samples <- unique(data$sample)
target_samples <- all_samples[all_samples != "MIgG"] # 排除 IgG
sample_pairs <- combn(target_samples, 2, simplify = FALSE)


for (pair in sample_pairs) {
  smp1 <- pair[1]
  smp2 <- pair[2]
  pair_name <- paste0(smp1, "_vs_", smp2)
  message(paste("Processing pair:", pair_name))
  
  pair_data <- data %>%
    filter(sample %in% pair & geneID %in% keep_genes & Family == "protein_coding") %>%
    group_by(geneID) %>%
    # 筛选：两个样本都大于 threshold，且两个样本都有数据
    filter(all(lfc_over_igg > threshold) & n() == 2) %>% 
    summarise(mean_lfc = mean(lfc_over_igg), .groups = 'drop') %>%
    arrange(desc(mean_lfc))

  intersect_genes_list <- pair_data$mean_lfc
  names(intersect_genes_list) <- pair_data$geneID

  message(paste("  - Found", length(intersect_genes_list), "shared genes."))

  if (length(intersect_genes_list) < 10) {
    message("  - Skipping: Too few shared genes.")
    next
  }
  
  gse_res <- tryCatch({
    gseGO(geneList = intersect_genes_list,
          OrgDb = org.Hs.eg.db,
          keyType = "ENSEMBL",
          ont = "BP",
          pvalueCutoff = 0.05,
          verbose = FALSE)
  }, error = function(e) return(NULL))
  
  #sig_genes <- smp_data %>% filter(lfc_over_igg > 1) %>% pull(geneID)

  # 4. 执行 GO 富集 (ORA)
  ego_res <- enrichGO(gene = names(intersect_genes_list),
                      universe = keep_genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pvalueCutoff = 0.05)

  # 5. 保存结果与绘图
  pair_out <- file.path(output_dir, pair_name)
  if (!dir.exists(pair_out)) dir.create(pair_out)
  
  write.csv(as.data.frame(gse_res), file.path(pair_out, paste0(pair_name, "_GSEA_results.csv")))
  write.csv(as.data.frame(ego_res), file.path(pair_out, paste0(pair_name, "_shared_GO.csv")))
  
  
  if (nrow(gse_res) > 0) {
    p1 <- dotplot(gse_res, showCategory=15) + 
      ggtitle(paste("Shared GO:", smp1, "&", smp2, "(t >", threshold, ")"))
    ggsave(file.path(pair_out, paste0(pair_name, "_GSEA_dotplot.pdf")), p1, width = 8, height = 6)
  }
  
  if (nrow(ego_res) > 0) {
    p2 <- barplot(ego_res, showCategory=15) + 
      ggtitle(paste("Shared GO:", smp1, "&", smp2, "(t >", threshold, ")"))
    ggsave(file.path(pair_out, paste0(pair_name, "_barplot.pdf")), p2, width = 8, height = 6)
  }
}



message("all done! check ", output_dir, "for results.")
