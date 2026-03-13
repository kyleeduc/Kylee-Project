############################################################################
# Author: Kylee Duczyminski
# About: Collects all of the genes in the adult mouse brain hippocampus that are
# considered imprinted in at least one mouse. This means that the allele was at
# least 85% expressed maternally or paternally.
# Input: data/processed/HIP_RNA-AlSp_paternalASE.csv
# Output: data/processed/HIP_RNA-mouse_imprinted_genes.csv
############################################################################

# Load libraries
library(tidyverse)
library(dplyr)
library(readr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project/")

# Read in data
hip_rna_ase <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE.csv")

# If a row contains 12 NA values, then delete the full row
hip_rna_ase <- hip_rna_ase[rowSums(is.na(hip_rna_ase)) != 12, ]

# If a row contains 12 "Not Imprinted" strings, then delete the full row
hip_rna_ase <- hip_rna_ase[rowSums(hip_rna_ase == "N/A", na.rm = TRUE) != 12, ]

# Create a new column that categorizes imprinting
hip_rna_ase <- hip_rna_ase %>%
  mutate(
    Imprinting_Category = case_when(
      # Check all imprinting_status columns for M, P, or N
      
      # Both M and P found
      (rowSums(across(ends_with("_ImprintStatus"), ~. == "Maternally Imprinted"), na.rm = TRUE) > 0) &
        (rowSums(across(ends_with("_ImprintStatus"), ~. == "Paternally Imprinted"), na.rm = TRUE) > 0) ~ "M&P",
      
      # Only M found
      (rowSums(across(ends_with("_ImprintStatus"), ~. == "Maternally Imprinted"), na.rm = TRUE) > 0) &
        (rowSums(across(ends_with("_ImprintStatus"), ~. == "Paternally Imprinted"), na.rm = TRUE) == 0) ~ "M",
      
      # Only P found
      (rowSums(across(ends_with("_ImprintStatus"), ~. == "Paternally Imprinted"), na.rm = TRUE) > 0) &
        (rowSums(across(ends_with("_ImprintStatus"), ~. == "Maternally Imprinted"), na.rm = TRUE) == 0) ~ "P",
      
      # Not imprinted (all are "Biallelic" or "N/A")
      (rowSums(across(ends_with("_ImprintStatus"), ~. == "Biallelic" | . == "N/A"), na.rm = TRUE) > 0) &
        (rowSums(across(ends_with("_ImprintStatus"), ~. != "Biallelic" & . != "N/A"), na.rm = TRUE) == 0) ~ "N",
      
      .default = "Unknown"
    )
  )

# View the results
print(hip_rna_ase %>% select(Gene_Name, Imprinting_Category) %>% head(20))

# Save the results to a new CSV file
hip_rna_ase %>%
  select(Gene_Name, Imprinting_Category) %>%
  write_csv("data/processed/HIP_RNA-mouse_imprinted_genes.csv")
