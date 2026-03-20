###############################################################################
# Author: Kylee Duczyminski
# About: This script is used to generate a heatmap using information about TPM for each gene
# Input: data/raw/HIP_RNA-TPM.xlsx
# Output: data/processed/HIP_RNA-TPM_heatmap.pdf
###############################################################################

# Load libraries
library(readxl)
library(dplyr)
library(readr)
library(pheatmap)
library(ggplot2)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Read your data
raw_data_TPM <- read_excel("data/raw/HIP_RNA-TPM.xlsx")

# Select sample columns
old_names <- c("1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250", 
               "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256")

new_names <- c("WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2", 
               "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2")

# Rename columns manually
colnames(raw_data_TPM)[match(old_names, colnames(raw_data_TPM))] <- new_names

# Organize the columns by genotype and number: WT_M, WT_F, HT_F, KO_M
ordered_columns <- c("WT_M1", "WT_M2", "WT_M3", "WT_F1", "WT_F2", "WT_F3", 
                     "HT_F1", "HT_F2", "HT_F3", "KO_M1", "KO_M2", "KO_M3")
raw_data_TPM <- raw_data_TPM %>% select(Gene_Name, all_of(ordered_columns))

# Select only the sample columns (exclude Gene_Name)
sample_columns <- raw_data_TPM %>% select(all_of(ordered_columns))

# Convert to matrix and set row names to gene names
data_matrix <- as.matrix(sample_columns)
rownames(data_matrix) <- raw_data_TPM$Gene_Name

# Replace NA with 0
data_matrix[is.na(data_matrix)] <- 0

# Log transform TPM
log_matrix <- log2(data_matrix + 1)

# Row-wise z-score normalization
data_matrix_normalized <- t(scale(t(log_matrix)))

# Replace any NA values caused by zero-variance rows
data_matrix_normalized[is.na(data_matrix_normalized)] <- 0

# Sort rows by WT_M1 expression (highest to lowest)
sort_order <- order(data_matrix_normalized[, "WT_M1"], decreasing = TRUE)
data_matrix_sorted <- data_matrix_normalized[sort_order, ]

# Save heatmap as PDF
# pdf("figures/HIP_RNA-TPM_heatmap.pdf", width = 10, height = 12)
# 
# heatmap(data_matrix_sorted,
#         main = "All Genes Expression Heatmap",
#         # Rowv = NA,
#         # Colv = NA,
#         scale = "none")
#         # margins = c(5, 10))
# 
# dev.off()

# # Normalize using z-score
# data_matrix_normalized <- scale(data_matrix)

# # Create heatmap using pheatmap
# pheatmap(data_matrix_normalized,
#          main = "All 26,000 Genes Expression Heatmap",
#          show_rownames = FALSE,
#          clustering_method = "ward.D2",
#          color = colorRampPalette(c("blue", "white", "red"))(50))

###############################################################################

# # Create a heatmap with heatmap()
# heatmap_generated <- heatmap(data_matrix_normalized,
#         main = "All 26,000 Genes Expression Heatmap",
#         Rowv = NA, 
#         Colv = NA, 
#         scale = "row",
#         margins = c(5, 10))
# 
# # Save heatmap_generated as a PDF to figures/
# pdf("figures/HIP_RNA-TPM_heatmap.pdf", width = 10, height = 10)
# heatmap(data_matrix_normalized,
#         main = "All 26,000 Genes Expression Heatmap",
#         Rowv = NA, 
#         Colv = NA, 
#         scale = "row",
#         margins = c(5, 10))
# dev.off()

gene_variance <- apply(log_matrix, 1, var)
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:500]

top_matrix <- log_matrix[top_genes, ]
top_matrix_scaled <- t(scale(t(top_matrix)))
top_matrix_scaled[is.na(top_matrix_scaled)] <- 0

# pheatmap(top_matrix_scaled,
#          clustering_method = "average",
#          show_rownames = FALSE)


# FURTHER analysis
sample_names <- colnames(data_matrix_normalized)

# Creates a data frame for column annotations based on sample Genotype and Sex
annotation_col <- data.frame(
  Genotype = ifelse(grepl("WT", sample_names), "WT",
                    ifelse(grepl("HT", sample_names), "HT", "KO")),
  Sex = ifelse(grepl("_M", sample_names), "Male", "Female")
)

# Set row names of annotation_col to sample names for proper alignment in pheatmap
rownames(annotation_col) <- sample_names

# Define colors for annotations
annotation_colors <- list(
  Genotype = c(WT = "#1f77b4", HT = "#ff7f0e", KO = "#2ca02c"),
  Sex = c(Male = "#4daf4a", Female = "#e41a1c")
)

# Generate and save heatmap with annotations
pdf("figures/HIP_RNA-TPM_heatmap.pdf", width = 10, height = 12)
pheatmap(top_matrix_scaled,
         clustering_method = "average",
         show_rownames = FALSE,
         cluster_rows = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors)
dev.off()

# PCA
# Find genes with non-zero variance across samples
gene_var <- apply(log_matrix, 1, var, na.rm = TRUE)
log_matrix_pca <- log_matrix[gene_var > 0, ]

# Run PCA
pca <- prcomp(t(log_matrix_pca), scale. = TRUE)

# Make PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Sample = colnames(log_matrix_pca),
  Genotype = annotation_col$Genotype,
  Sex = annotation_col$Sex
)

library(ggplot2)

# Generate and save PCA plot as a .pdf
pdf("figures/HIP_RNA-TPM_PCA.pdf", width = 8, height = 8)
ggplot(pca_df, aes(PC1, PC2, color = Genotype, shape = Sex)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1) +
  theme_minimal() +
  ggtitle("PCA of RNA-seq Samples")
dev.off()

