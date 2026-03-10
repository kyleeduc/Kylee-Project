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
  rename(Expressed_Allele = Expressed.Allele)  # Rename to match imprinted_genes

imprinted_genes <- read.csv("data/raw/imprinted_gene_list.csv") %>%
  mutate(source = "imprinted_gene_list")

#Combine the datasets
imprinted_reference <- bind_rows(imprinted_brain, imprinted_genes) %>%
  group_by(Gene) %>%
  summarise(
    source = paste(unique(source), collapse = " & "),
    Expressed_Allele = first(Expressed_Allele, na_rm = TRUE),
    Status = first(Status, na_rm = TRUE),
    .groups = "drop"
  ) %>%
  distinct(Gene, .keep_all = TRUE)

# Rename sample name columns for easier handling using regex with word boundaries
# WT F are identified as samples X1, X8, and X9
# WT M are identified as samples X5, X6, and X11
# Instead of sequential gsub calls, use a mapping approach
# Map sample IDs to their new names
sample_mapping <- c(
  "X1" = "WT_F1",
  "X8" = "WT_F2",
  "X9" = "WT_F3",
  "X5" = "WT_M1",
  "X6" = "WT_M2",
  "X11" = "WT_M3"
)

# Rename columns by replacing the sample ID prefix
colnames(hip_rna_ase) <- sapply(colnames(hip_rna_ase), function(col) {
  for (old_id in names(sample_mapping)) {
    if (startsWith(col, paste0(old_id, "_"))) {
      return(sub(paste0("^", old_id, "_"), paste0(sample_mapping[old_id], "_"), col))
    }
  }
  return(col)  # Return unchanged if no match
})

# Verify the renaming worked
colnames(hip_rna_ase)

# Delete any rows from the hip_rna_ase that contain NA values
hip_rna_ase <- hip_rna_ase %>%
  drop_na()

# Combine imprinted_brain and imprinted_genes, keeping imprinted_brain for duplicates
imprinted_reference <- bind_rows(imprinted_brain, imprinted_genes) %>%
  distinct(Gene, .keep_all = TRUE)  # Keep first occurrence (imprinted_brain) if duplicate

# Filter the imprinted reference data to only include genes that are in our ASE dataset
imprinted_reference_filtered <- imprinted_reference %>%
  filter(Gene %in% hip_rna_ase$base_gene_name)

# Merge the imprinted reference data with our ASE data
imprinting_analysis <- hip_rna_ase %>%
  inner_join(imprinted_reference_filtered, by = c("base_gene_name" = "Gene"))

# Calculate averages and consistency
imprinting_analysis <- imprinting_analysis %>%
  rowwise() %>%
  mutate(
    # Trim whitespace from Expressed.Allele
    Expressed_Allele = trimws(Expressed_Allele),
    WT_F_avg = mean(c(WT_F1_S245_allelic_expression_ratio, 
                      WT_F2_S252_allelic_expression_ratio, 
                      WT_F3_S253_allelic_expression_ratio), na.rm = TRUE),
    WT_M_avg = mean(c(WT_M1_S249_allelic_expression_ratio, 
                      WT_M2_S250_allelic_expression_ratio, 
                      WT_M3_S255_allelic_expression_ratio), na.rm = TRUE),
    WT_Overall_avg = mean(c(WT_F_avg, WT_M_avg), na.rm = TRUE),
    Expression_Pattern = case_when(
      Expressed_Allele == "Maternal" & WT_Overall_avg >= 0.7 ~ "Matches Maternal",
      Expressed_Allele == "Paternal" & WT_Overall_avg <= -0.7 ~ "Matches Paternal",
      TRUE ~ "Does Not Match"
    ),
    # Check if all WT females are in the same range
    WT_Female_Consistent = case_when(
      (WT_F1_S245_allelic_expression_ratio >= 0.7 & 
         WT_F2_S252_allelic_expression_ratio >= 0.7 & 
         WT_F3_S253_allelic_expression_ratio >= 0.7) |
        (WT_F1_S245_allelic_expression_ratio > -0.7 & WT_F1_S245_allelic_expression_ratio < 0.7 &
           WT_F2_S252_allelic_expression_ratio > -0.7 & WT_F2_S252_allelic_expression_ratio < 0.7 &
           WT_F3_S253_allelic_expression_ratio > -0.7 & WT_F3_S253_allelic_expression_ratio < 0.7) |
        (WT_F1_S245_allelic_expression_ratio <= -0.7 & 
           WT_F2_S252_allelic_expression_ratio <= -0.7 & 
           WT_F3_S253_allelic_expression_ratio <= -0.7) ~ "Yes",
      TRUE ~ "No"
    ),
    # Check if all WT males are in the same range
    WT_Male_Consistent = case_when(
      (WT_M1_S249_allelic_expression_ratio >= 0.7 & 
         WT_M2_S250_allelic_expression_ratio >= 0.7 & 
         WT_M3_S255_allelic_expression_ratio >= 0.7) |
        (WT_M1_S249_allelic_expression_ratio > -0.7 & WT_M1_S249_allelic_expression_ratio < 0.7 &
           WT_M2_S250_allelic_expression_ratio > -0.7 & WT_M2_S250_allelic_expression_ratio < 0.7 &
           WT_M3_S255_allelic_expression_ratio > -0.7 & WT_M3_S255_allelic_expression_ratio < 0.7) |
        (WT_M1_S249_allelic_expression_ratio <= -0.7 & 
           WT_M2_S250_allelic_expression_ratio <= -0.7 & 
           WT_M3_S255_allelic_expression_ratio <= -0.7) ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  ungroup()

write.csv(imprinting_analysis, "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv", row.names = FALSE)

