#!/usr/bin/env Rscript

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

# Perform PCA (exclude 'TE' column as it's a label)
expr_data_filtered <- expr_data[apply(expr_data, 1, var) > 0, ]
pca_result <- prcomp(t(expr_data_filtered), scale. = TRUE)

# Extract the variance explained by each principal component
explained_variance <- summary(pca_result)$importance[2, ]  # Get the proportion of variance

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




# Perform hierarchical clustering
dist_matrix <- dist(t(expr_data))  # Calculate the distance matrix (transpose to get samples as rows)
hc_result <- hclust(dist_matrix, method = "complete")  # Perform hierarchical clustering using 'complete' linkage

# Create dendrogram data
dendro_data <- dendro_data(hc_result)

# Extract leaf node labels and their coordinates (x, y positions)
leaf_nodes <- dendro_data$labels$label
leaf_positions <- dendro_data$labels[, c("x", "y")]

# Check if the number of leaf positions matches the number of leaf nodes
if (length(leaf_nodes) != nrow(leaf_positions)) {
  stop("The number of leaf nodes and leaf positions do not match!")
}

# Create a data frame for leaf labels
leaf_labels <- data.frame(x = leaf_positions$x, y = leaf_positions$y, label = leaf_nodes)

# Plot the hierarchical clustering dendrogram with sample labels
hc_plot <- ggplot() + 
  geom_segment(data = segment(dendro_data), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = leaf_labels, aes(x = x, y = y, label = label), vjust = 1.5, hjust = -0.1, size = 3, angle=45) +  # Add sample names
  theme_minimal() +
  ggtitle("Hierarchical Clustering of TE Expression")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save the hierarchical clustering plot with sample names
ggsave(file.path(file_dir, "hc_plot.png"), plot = hc_plot)

# Optionally, save the hierarchical clustering tree as a CSV for further inspection
# write.csv(hc_result$merge, file.path(file_dir, "hc_results.csv"))


# --- Hierarchical Clustering on PCA-reduced data ---
# Select the first two principal components for clustering (or you can choose more)
#pca_reduced_data <- pca_result$x[, 1:2]  # Use first two principal components
#pca_reduced_data <- pca_result$x

explained_variance <- summary(pca_result)$importance[2, ]  # Variance explained by each PC
cumulative_variance <- cumsum(explained_variance)  # Cumulative variance explained
# Find the number of PCs that explain more than 80% of the variance
n_pcs <- which(cumulative_variance >= 0.80)[1]  # Find the first PC where the cumulative variance exceeds 80%
# Select the PCs that explain 80% of the variance
pca_reduced_data <- pca_result$x[, 1:n_pcs] 

cat("Number of PCs that explain more than 80% of the variance:", n_pcs, "\n")
cat("Variance explained by these", n_pcs, "PCs:", cumulative_variance[n_pcs], "\n")

# Perform hierarchical clustering on the PCA-reduced data
dist_matrix_pca <- dist(pca_reduced_data)  # Calculate the distance matrix on PCA-reduced data
hc_result_pca <- hclust(dist_matrix_pca, method = "complete")  # Perform hierarchical clustering using 'complete' linkage

# Create dendrogram data for PCA-based HC
dendro_data_pca <- dendro_data(hc_result_pca)

# Extract leaf node labels and their coordinates (x, y positions) for PCA-based HC
leaf_nodes_pca <- dendro_data_pca$labels$label
leaf_positions_pca <- dendro_data_pca$labels[, c("x", "y")]

# Check if the number of leaf positions matches the number of leaf nodes
if (length(leaf_nodes_pca) != nrow(leaf_positions_pca)) {
  stop("The number of leaf nodes and leaf positions for PCA data do not match!")
}

# Create a data frame for leaf labels for PCA-based HC
leaf_labels_pca <- data.frame(x = leaf_positions_pca$x, y = leaf_positions_pca$y, label = leaf_nodes_pca)

# Plot the hierarchical clustering dendrogram on PCA-reduced data with sample labels
hc_plot_pca <- ggplot() + 
  geom_segment(data = segment(dendro_data_pca), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = leaf_labels_pca, aes(x = x, y = y, label = label), vjust = 1.5, hjust = -0.1, size = 3,angle=45) +  # Add sample names
  theme_minimal() +
  ggtitle("Hierarchical Clustering of PCA-Reduced TE Expression Data")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the hierarchical clustering plot for PCA-reduced data
ggsave(file.path(file_dir, "hc_plot_pca.png"), plot = hc_plot_pca)

# Optionally, save the PCA-based hierarchical clustering tree as a CSV for further inspection
# write.csv(hc_result_pca$merge, file.path(file_dir, "hc_results_pca.csv"))
