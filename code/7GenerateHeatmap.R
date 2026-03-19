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

# Normalize using z-score
data_matrix_normalized <- scale(data_matrix)

# # Create heatmap using pheatmap
# pheatmap(data_matrix_normalized,
#          main = "All 26,000 Genes Expression Heatmap",
#          show_rownames = FALSE,
#          clustering_method = "ward.D2",
#          color = colorRampPalette(c("blue", "white", "red"))(50))

###############################################################################

# Create a heatmap with heatmap()
heatmap_generated <- heatmap(data_matrix_normalized,
        main = "All 26,000 Genes Expression Heatmap",
        Rowv = NA, 
        Colv = NA, 
        scale = "row",
        margins = c(5, 10))

# Save the heatmap as a PDF
pdf("data/processed/HIP_RNA-TPM_heatmap.pdf", width = 8, height = 12)
