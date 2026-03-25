###############################################################################
# Author: Kylee Duczyminski
# About: This script will generate a heatmap for the 112 published imprinted genes
# in adult mice.
# Input: data/raw/HIP_RNA-TPM.xlsx
#        data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
# Output: figures/HIP_RNA-TPM_heatmap_112_published_genes.pdf
###############################################################################

# Load libraries
library(tidyverse)
library(readxl)
library(pheatmap)
library(readr)
library(dplyr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load data
tpm <- read_excel("data/raw/HIP_RNA-TPM.xlsx")
imprinted_genes <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# If Gene_Name from imprinted_genes is present in tpm, then keep that row in tpm
tpm_imprinted <- tpm %>%
  filter(Gene_Name %in% imprinted_genes$Gene_Name)

# Select sample columns
old_names <- c("1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250", 
               "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256")

new_names <- c("WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2", 
               "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2")

# Rename columns manually
colnames(tpm_imprinted)[match(old_names, colnames(tpm_imprinted))] <- new_names

# Organize the columns by genotype and number: WT_M, WT_F, HT_F, KO_M
ordered_columns <- c("WT_M1", "WT_M2", "WT_M3", "WT_F1", "WT_F2", "WT_F3", 
                     "HT_F1", "HT_F2", "HT_F3", "KO_M1", "KO_M2", "KO_M3")
tpm_imprinted <- tpm_imprinted %>% select(Gene_Name, all_of(ordered_columns))

# Prepare data for heatmap
heatmap_data <- tpm_imprinted %>%
  column_to_rownames("Gene_Name") %>%
  as.matrix()

# # Replace any NA, NaN, or Inf values with 0
# heatmap_data[is.na(heatmap_data)] <- 0
# heatmap_data[is.infinite(heatmap_data)] <- 0

# Replace any 0s in tpm with a value of NA
heatmap_data[heatmap_data == 0] <- NA

# Log-transform TPM values
heatmap_data <- log2(heatmap_data)

# # Filter out genes with all zeros or mostly zeros
# heatmap_data <- heatmap_data[rowSums(heatmap_data) > 0, ]

# Change NA values to 0
heatmap_data[is.na(heatmap_data)] <- 0

# Remove rows with zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, var) != 0, ]

# Create heatmap
pdf("figures/published_imprinted_genes_heatmap.pdf", width = 12, height = 20)
pheatmap(heatmap_data,
         main = "Imprinted Genes Expression Heatmap",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         cellwidth = 30,
         cellheight = 5,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = "row")
dev.off()

# Heatmap that keeps the genes with zero expression at the bottom

# # Separate genes into two groups
# heatmap_data_nonzero <- heatmap_data[rowSums(heatmap_data) > 0, ]
# heatmap_data_zero <- heatmap_data[rowSums(heatmap_data) == 0, ]
# 
# # Combine them back
# heatmap_data_combined <- rbind(heatmap_data_nonzero, heatmap_data_zero)
# 
# # Create heatmap
# pheatmap(heatmap_data_combined,
#          main = "Imprinted Genes Expression Heatmap",
#          cluster_rows = FALSE,  # Disable clustering to keep groups separate
#          cluster_cols = FALSE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize = 6,
#          fontsize_row = 5,
#          cellwidth = 30,
#          cellheight = 2,
#          color = colorRampPalette(c("blue", "white", "red"))(50),
         # scale = "row")
