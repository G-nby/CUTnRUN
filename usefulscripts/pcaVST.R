#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  #library(factoextra)
  library(cluster)
})

library(ggplot2)  
library(tools)
library(ggdendro)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Usage: Rscript pcatmp.R <file/path/EXPR.csv>")
}

# Read the expression data
data_file <- args[1]
expr_data <- read.csv(data_file, row.names = 1)

# Check if the data is loaded correctly
if (is.null(expr_data)) {
  stop("Failed to load expression data.")
}

coldata <- data.frame(sample = colnames(expr_data), row.names = colnames(expr_data))
dds <- DESeqDataSetFromMatrix(countData = expr_data,
                              colData = coldata,
                              design = ~ 1)
                              
dds <- dds[rowSums(counts(dds)) > 10, ]

# blind vst, for no group preset
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

pca_result <- prcomp(t(vsd_mat))
explained_variance <- summary(pca_result)$importance[2, ]


# Perform PCA (exclude 'TE' column as it's a label)
#expr_data_filtered <- expr_data[apply(expr_data, 1, var) > 0, ]
#pca_result <- prcomp(t(expr_data_filtered), scale. = TRUE)

# Extract the variance explained by each principal component
#explained_variance <- summary(pca_result)$importance[2, ]  # Get the proportion of variance

pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Sample = colnames(expr_data))

# Add variance explained to axis labels
x_label <- paste("PC1 (", round(explained_variance[1] * 100, 2), "%)", sep = "")
y_label <- paste("PC2 (", round(explained_variance[2] * 100, 2), "%)", sep = "")

# Print the total variance explained by PCA
total_variance_explained <- sum(explained_variance) * 100  # Total percentage variance explained by all PCs
cat("Total variance explained by PCA: ", round(total_variance_explained, 2), "%\n")
  
# Plot the PCA results
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample, label = Sample)) +
  geom_point(size = 3) +  # Size of points
  geom_text(vjust = -1) +  # Add labels above the points
  ggtitle("PCA of TE Expression") +
  xlab(x_label) +
  ylab(y_label) +
  theme_minimal()

file_dir <- dirname(data_file)
ggsave(file.path(file_dir, "pca_plot.png"), plot = pca_plot)

# Save PCA numerical results to a CSV file
pca_values <- data.frame(Sample = colnames(expr_data), pca_result$x)  # Add sample names to the PCA result
write.csv(pca_values, file.path(file_dir, "pca_results.csv"), row.names = FALSE)




