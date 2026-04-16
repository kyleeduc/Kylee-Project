##############################################################################
# Author: Kylee Duczyminski
# About: This script will create a heatmap of genes that are considered imprinted, 
# and identify samples  that have slight biases toward one parental allele. If
# samples have no bias, they will be colored white, if maternal bias = red, if 
# paternal bias = blue. If NA value for the sample exists, color the sample gray.
# Input: HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv from data/processed
# Output: heatmap of imprinted genes and samples with bias toward one parental allele
##############################################################################

# Load libraries
library(tidyverse)
library(readr)
library(dplyr)
library(pheatmap)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Read in data
imprinting_data <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# Organize samples by genotype
ordered_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3",
  "HT_F1", "HT_F2", "HT_F3",
  "KO_M1", "KO_M2", "KO_M3"
)

# Build the matching AER column names
ordered_aer_columns <- paste0(ordered_samples, "_AER")

# Keep Gene_Name + reordered AER columns
imprinting_data <- imprinting_data %>%
  select(Gene_Name, all_of(ordered_aer_columns))

# Convert to matrix (genes as rows)
mat <- imprinting_data %>%
  column_to_rownames("Gene_Name") %>%
  as.matrix()

# ------------------------------------------------------------------
# 1. Categorize AER values
# ------------------------------------------------------------------
categorize_ratio <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x >= 0.7 ~ "red",
    x <= -0.7 ~ "blue",
    x > 0 & x < 0.7 ~ "light_red",
    x > -0.7 & x < 0 ~ "light_blue",
    x == 0 ~ "white"
  )
}

cat_mat <- matrix(
  categorize_ratio(as.vector(mat)),
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
)

# ------------------------------------------------------------------
# 2. Convert categories → numbers (required for pheatmap)
# ------------------------------------------------------------------
category_to_number <- c(
  "blue" = 1,
  "light_blue" = 2,
  "white" = 3,
  "light_red" = 4,
  "red" = 5,
  "NA" = 6
)

plot_mat <- matrix(
  category_to_number[cat_mat],
  nrow = nrow(cat_mat),
  ncol = ncol(cat_mat),
  dimnames = dimnames(cat_mat)
)

# ------------------------------------------------------------------
# 3. Define exact colors
# ------------------------------------------------------------------
colors <- c(
  "blue",
  "lightskyblue",
  "white",
  "lightcoral",
  "red",
  "gray"
)

breaks <- seq(0.5, 6.5, by = 1)

# ------------------------------------------------------------------
# 4. Add genotype annotation (optional but VERY helpful visually)
# ------------------------------------------------------------------
annotation_col <- data.frame(
  Genotype = c(
    rep("WT_M", 3),
    rep("WT_F", 3),
    rep("HT_F", 3),
    rep("KO_M", 3)
  )
)
rownames(annotation_col) <- colnames(plot_mat)

# Optional: color the genotype labels
annotation_colors <- list(
  Genotype = c(
    WT_M = "#1f77b4",
    WT_F = "#ff7f0e",
    HT_F = "#2ca02c",
    KO_M = "#d62728"
  )
)

# ------------------------------------------------------------------
# 5. Plot and save heatmap
# ------------------------------------------------------------------

pheatmap(
  plot_mat,
  color = colors,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = "black",
  filename = "figures/heatmap_112_genes_tall.pdf",
  width = 8,
  height = 18
)

# Run a loop that calculates the coefficient of variation (SD/mean)
# of _AER values for each gene across samples within each genotype
genotypes <- c("WT_M", "WT_F", "HT_F", "KO_M")
cv_results <- data.frame(Gene_Name = rownames(mat))

for (genotype in genotypes) {
  aer_cols <- grep(paste0("^", genotype, "[0-9]+_AER$"), colnames(mat), value = TRUE)
  
  print(genotype)
  print(aer_cols)
  
  cv_results[[paste0(genotype, "_CV")]] <- apply(mat[, aer_cols, drop = FALSE], 1, function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    m <- mean(x, na.rm = TRUE)
    if (m == 0) {
      return(NA)
    }
    sd(x, na.rm = TRUE) / m
  })
}



