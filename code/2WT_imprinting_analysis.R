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

# Count the number of mice that are either "Maternally Imprinted" or "Paternally Imprinted" under columns ending in "_ImprintStatus" and create a column for this value
imprinting_analysis <- imprinting_analysis %>%
  rowwise() %>%
  mutate(Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) %in% c("Maternally Imprinted", "Paternally Imprinted"), na.rm = TRUE)) %>%
  ungroup()

# Count the number of mice that have "Maternally Imprinted" under the columns with the ending "_ImprintStatus" and create a column for this value
imprinting_analysis <- imprinting_analysis %>%
  rowwise() %>%
  mutate(Maternally_Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) == "Maternally Imprinted", na.rm = TRUE)) %>%
  ungroup()

# Count the number of mice that have "Paternally Imprinted" under the columns with the ending "_ImprintStatus" and create a column for this value
imprinting_analysis <- imprinting_analysis %>%
  rowwise() %>%
  mutate(Paternally_Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) == "Paternally Imprinted", na.rm = TRUE)) %>%
  ungroup()

# If Maternally_Imprinted_Mice_Count is equal to Imprinted_Mice_Count, then create a column called "Overall_Imprinting_Status" and set it to "Maternally Imprinted". If Paternally_Imprinted_Mice_Count is equal to Imprinted_Mice_Count, then set "Overall_Imprinting_Status" to "Paternally Imprinted". If neither of these conditions are met, set "Overall_Imprinting_Status" to "Inconsistent Imprinting".
imprinting_analysis <- imprinting_analysis %>%
  mutate(Overall_Expression_Status = case_when(
    Maternally_Imprinted_Mice_Count == Imprinted_Mice_Count & Imprinted_Mice_Count > 0 ~ "Paternal",
    Paternally_Imprinted_Mice_Count == Imprinted_Mice_Count & Imprinted_Mice_Count > 0 ~ "Maternal",
    Imprinted_Mice_Count == 0 ~ "Not Imprinted",
    TRUE ~ "Inconsistent"
  ))

# If Overall_Expression_Status matches "Expressed_Allele", then create a column called "Data_Matching" and set it to "Yes". If Overall_Expression_Status does not match "Expressed_Allele", then set "Data_Matching" to "No". If Overall_Expression_Status is "Inconsistent", then set "Data_Matching" to "Inconsistent". 
imprinting_analysis <- imprinting_analysis %>%
  mutate(Data_Matching = case_when(
    Overall_Expression_Status == Expressed_Allele ~ "Yes",
    Overall_Expression_Status == "Inconsistent" ~ "Inconsistent",
    TRUE ~ "No"
  ))

# Calculate the percent at which the Overall_Expression_Status matches the Expressed_Allele for the genes that are not "Inconsistent". If they are "Inconsistent," check "Expressed_Allele" to identify which expression pattern is expected. If "Maternal" is expected, use the formula (Paternally_Imprinted_Mice_Count/Imprinted_Mice_Count) * 100. If "Paternal" is expected, use the formula (Maternally_Imprinted_Mice_Count/Imprinted_Mice_Count) * 100. Create a column called "Percent_Matching" and set it to this value.
imprinting_analysis <- imprinting_analysis %>%
  mutate(Expressed_Allele = trimws(Expressed_Allele)
    ,Percent_Matching = case_when(
    Data_Matching == "Yes" ~ 100,
    Data_Matching == "No" ~ 0,
    Data_Matching == "Inconsistent" & Expressed_Allele == "Maternal" ~ round((Paternally_Imprinted_Mice_Count / Imprinted_Mice_Count) * 100, 1),
    Data_Matching == "Inconsistent" & Expressed_Allele == "Paternal" ~ round((Maternally_Imprinted_Mice_Count / Imprinted_Mice_Count) * 100, 1),
    TRUE ~ NA_real_
  ))

# Create a column titled Percent_Mismatch and set it to 100 - Percent_Matching
imprinting_analysis <- imprinting_analysis %>%
  mutate(Percent_Mismatch = case_when(
    Overall_Expression_Status == "Not Imprinted" ~ 0,
    TRUE ~ round(100 - Percent_Matching, 1)
  ))

write.csv(imprinting_analysis, "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv", row.names = FALSE)

