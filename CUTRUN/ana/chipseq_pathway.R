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
  stop("用法: Rscript script.R <input_dir> [output_dir]")
}

input_dir <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], file.path(input_dir, "pathway_analysis"))

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

#print(target_samples)

for (smp in target_samples) {
  message(paste("Processing sample:", smp))
  
  #print()
  smp_data <- data %>% 
    filter(sample == smp & geneID %in% keep_genes & Family == "protein_coding") %>%
    arrange(desc(lfc_over_igg))
  print(smp_data)
  
  if (nrow(smp_data) == 0) next
  
  gene_list <- smp_data$lfc_over_igg
  names(gene_list) <- smp_data$geneID
  
  head(gene_list)
  
  gse_res <- gseGO(geneList = gene_list,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL",
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   verbose = FALSE)
  
  sig_genes <- smp_data %>% filter(lfc_over_igg > 1) %>% pull(geneID)
  
  ego_res <- enrichGO(gene = sig_genes,
                      universe = keep_genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pvalueCutoff = 0.05)
  
  sample_out <- file.path(output_dir, smp)
  if (!dir.exists(sample_out)) dir.create(sample_out)
  
  write.csv(as.data.frame(gse_res), file.path(sample_out, paste0(smp, "_GSEA_results.csv")))
  write.csv(as.data.frame(ego_res), file.path(sample_out, paste0(smp, "_GO_ORA_results.csv")))
  
  if (nrow(gse_res) > 0) {
    p1 <- dotplot(gse_res, showCategory=15) + ggtitle(paste(smp, "GSEA BP"))
    ggsave(file.path(sample_out, paste0(smp, "_GSEA_dotplot.pdf")), p1, width = 8, height = 6)
  }
  
  if (nrow(ego_res) > 0) {
    p2 <- barplot(ego_res, showCategory=15) + ggtitle(paste(smp, "GO Enrichment"))
    ggsave(file.path(sample_out, paste0(smp, "_GO_barplot.pdf")), p2, width = 8, height = 6)
  }
  
  message(paste("Sample", smp, "done!"))
}

message("all done! check ", output_dir, "for results.")
