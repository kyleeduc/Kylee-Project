###############################################################################
# Author: Kylee Duczyminski
# About: This script performs an analysis of our WT mice to determine if they
# follow the expected imprinting patterns.
# Input: data/processed/HIP_RNA-AlSp_parentalASE.csv
#        data/raw/imprinted_brain_jf1x129s1.csv
#        data/raw/imprinted_gene_list.csv
# Output: data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
###############################################################################

# Set proper working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)

# Load in the data
hip_rna_ase <- read.csv("data/processed/HIP_RNA-AlSp_parentalASE.csv")

imprinted_brain <- read.csv("data/raw/imprinted_brain_jf1x129s1.csv") %>%
  mutate(source = "imprinted_brain") %>%
  select(-Strand) %>%
  rename(Expressed_Allele = Expressed.Allele)

imprinted_genes <- read.csv("data/raw/imprinted_gene_list.csv") %>%
  mutate(source = "imprinted_gene_list")

# Keep rows that have at least one non-NA AER value
hip_rna_ase <- hip_rna_ase %>%
  filter(if_any(ends_with("_AER"), ~!is.na(.)))

nrow(hip_rna_ase)

# Combine imprinted_brain and imprinted_genes, keeping imprinted_brain for duplicates
imprinted_reference <- bind_rows(imprinted_brain, imprinted_genes) %>%
  distinct(Gene, .keep_all = TRUE)  # Keep first occurrence (imprinted_brain) if duplicate

# Filter the imprinted reference data to only include genes that are in our ASE dataset
imprinted_reference_filtered <- imprinted_reference %>%
  filter(Gene %in% hip_rna_ase$Gene_Name)

# Merge the imprinted reference data with our ASE data
imprinting_analysis <- hip_rna_ase %>%
  inner_join(imprinted_reference_filtered, by = c("Gene_Name" = "Gene"))

# Calculate averages and consistency
imprinting_analysis <- imprinting_analysis %>%
  rowwise() %>%
  mutate(
    # Trim whitespace from Expressed.Allele
    Expressed_Allele = trimws(Expressed_Allele),
    WT_F_avg = mean(c(WT_F1_AER, 
                      WT_F2_AER, 
                      WT_F3_AER), na.rm = TRUE),
    WT_M_avg = mean(c(WT_M1_AER, 
                      WT_M2_AER, 
                      WT_M3_AER), na.rm = TRUE),
    WT_Overall_avg = mean(c(WT_F_avg, WT_M_avg), na.rm = TRUE),
    Expression_Pattern = case_when(
      Expressed_Allele == "Maternal" & WT_Overall_avg >= 0.7 ~ "Matches Maternal",
      Expressed_Allele == "Paternal" & WT_Overall_avg <= -0.7 ~ "Matches Paternal",
      TRUE ~ "Does Not Match"
    ),
    # Check if all WT females are in the same range
    WT_Female_Consistent = case_when(
      (WT_F1_AER >= 0.7 & 
         WT_F2_AER >= 0.7 & 
         WT_F3_AER >= 0.7) |
        (WT_F1_AER > -0.7 & WT_F1_AER < 0.7 &
           WT_F2_AER > -0.7 & WT_F2_AER < 0.7 &
           WT_F3_AER > -0.7 & WT_F3_AER < 0.7) |
        (WT_F1_AER <= -0.7 & 
           WT_F2_AER <= -0.7 & 
           WT_F3_AER <= -0.7) ~ "Yes",
      TRUE ~ "No"
    ),
    # Check if all WT males are in the same range
    WT_Male_Consistent = case_when(
      (WT_M1_AER >= 0.7 & 
         WT_M2_AER >= 0.7 & 
         WT_M3_AER >= 0.7) |
        (WT_M1_AER > -0.7 & WT_M1_AER < 0.7 &
           WT_M2_AER > -0.7 & WT_M2_AER < 0.7 &
           WT_M3_AER > -0.7 & WT_M3_AER < 0.7) |
        (WT_M1_AER <= -0.7 & 
           WT_M2_AER <= -0.7 & 
           WT_M3_AER <= -0.7) ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  ungroup()

write.csv(imprinting_analysis, "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv", row.names = FALSE)

