###############################################################################
# Author: Kylee Duczyminski
# About: Calculate consistency between WT_F and WT_M and all WT by the final 
# imprinting status within their respective genotypes. Further analysis of HT
# and KO genotypes will be performed in this script.
# Input: data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
# Output: data/processed/genotype_consistency.csv
###############################################################################

# Load libraries
library(tidyverse)
library(readr)
library(dplyr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Read in data containing the 112 imprinted genes
imprinted_gene_data <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# Create a new data frame to save genotype-specific results. Columns include: "Gene_Name" & "Published_Status"
genotype_consistency <- imprinted_gene_data %>%
  select(Gene_Name, Expressed_Allele)

# Rename Expressed_Allele to Published_Expression
colnames(genotype_consistency)[2] <- "Published_Expression"

# Create a loop for each genotype
genotypes <- c("WT_F", "WT_M", "HT_F", "KO_M")

for (genotype in genotypes) {
  # Create a new column that reads in the number of samples that are either "Maternally Imprinted" or "Paternally Imprinted"
  imprinted_count_col <- paste0(genotype, "_Imprinted_Count")
  genotype_consistency[[imprinted_count_col]] <- imprinted_gene_data %>%
    select(starts_with(genotype)) %>%
    apply(1, function(x) sum(x == "Maternally Imprinted" | x == "Paternally Imprinted", na.rm = TRUE))
  
  # Create a column that counts "Maternally Imprinted" samples
  maternally_col <- paste0(genotype, "_Maternally_Imprinted_Count")
  genotype_consistency[[maternally_col]] <- imprinted_gene_data %>%
    select(starts_with(genotype)) %>%
    apply(1, function(x) sum(x == "Maternally Imprinted", na.rm = TRUE))
  
  # Create a column that counts "Paternally Imprinted" samples
  paternally_col <- paste0(genotype, "_Paternally_Imprinted_Count")
  genotype_consistency[[paternally_col]] <- imprinted_gene_data %>%
    select(starts_with(genotype)) %>%
    apply(1, function(x) sum(x == "Paternally Imprinted", na.rm = TRUE))
  
  # Create a column that states the Imprint_Status for this genotype
  imprint_status_col <- paste0(genotype, "_Imprint_Status")
  genotype_consistency[[imprint_status_col]] <- ifelse(
    genotype_consistency[[imprinted_count_col]] == 0, "Not Imprinted",
    ifelse(
      genotype_consistency[[maternally_col]] == genotype_consistency[[imprinted_count_col]], "Maternal",
      ifelse(
        genotype_consistency[[paternally_col]] == genotype_consistency[[imprinted_count_col]], "Paternal", "Inconsistent"
      )
    )
  )
  
  # Create a column for Percent_Matching
  percent_matching_col <- paste0(genotype, "_Percent_Matching")
  genotype_consistency <- genotype_consistency %>%
    mutate(
      Published_Expression = trimws(Published_Expression),
      !!sym(percent_matching_col) := case_when(
        !!sym(imprinted_count_col) == 0 ~ NA_real_,
        Published_Expression == "Maternal" ~ round((!!sym(paternally_col) / !!sym(imprinted_count_col)) * 100, 1),
        Published_Expression == "Paternal" ~ round((!!sym(maternally_col) / !!sym(imprinted_count_col)) * 100, 1),
        TRUE ~ NA_real_
      )
    )
}

# Hide the calculation columns, keeping only Gene_Name, Published_Expression, Imprint_Status columns, and Percent_Matching columns
genotype_consistency <- genotype_consistency %>%
  select(Gene_Name, Published_Expression, ends_with("_Imprint_Status"), ends_with("_Percent_Matching"))

# Write the final data frame to a new CSV file
write_csv(genotype_consistency, "data/processed/genotype_consistency.csv")







