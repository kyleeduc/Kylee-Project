###########################################################################
# Author: Kylee Duczyminski
# About: Cleaning up .csv file for easier future applications
# Input: data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
# Output: data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis_clean.csv
###########################################################################

# Set proper working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load necessary libraries
library(dplyr)
library(readr)

# Load in the data
imprinting_analysis <- read.csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# Rename columns ending in "_allelic_expression_ratio" to end in "_AER" for brevity
colnames(imprinting_analysis) <- gsub("_allelic_expression_ratio", "_AER", colnames(imprinting_analysis))
colnames(imprinting_analysis) <- gsub("_imprinting_status", "_ImprintStatus", colnames(imprinting_analysis))

# Rename columns to be more concise and list the mouse genotypes

# Rename columns starting with "X10_S254_" to "5cKO_M1"
colnames(imprinting_analysis) <- gsub("X10_S254_", "5cKO_M1_", colnames(imprinting_analysis))

# Rename columns starting with "X12_S256_" to "5cKO_M2"
colnames(imprinting_analysis) <- gsub("X12_S256_", "5cKO_M2_", colnames(imprinting_analysis))

# Rename columns starting with "X2_S246_" to "5cHT_F1"
colnames(imprinting_analysis) <- gsub("X2_S246_", "5cHT_F1_", colnames(imprinting_analysis))

# Rename columns starting with "X3_S247_" to "5cHT_F2"
colnames(imprinting_analysis) <- gsub("X3_S247_", "5cHT_F2_", colnames(imprinting_analysis))

# Rename columns starting with "X4_S248_" to "5cKO_M3"
colnames(imprinting_analysis) <- gsub("X4_S248_", "5cKO_M3_", colnames(imprinting_analysis))

# Rename columns starting with "X7_S251_" to "5cHT_F3"
colnames(imprinting_analysis) <- gsub("X7_S251_", "5cHT_F3_", colnames(imprinting_analysis))

# Print out column titles to check for correctness
print(colnames(imprinting_analysis))

# Delete all columns ending in "_AER"
imprinting_analysis_clean <- imprinting_analysis %>%
  select(-ends_with("_AER"))


# If value is "Not Imprinted" in columns ending in "_Imprinted_Status", change to "N"
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(across(ends_with("_ImprintStatus"), ~ ifelse(. == "Not Imprinted", "N", .)))

# If value is "Paternally Imprinted" in columns ending in "_Imprinted_Status", change to "P"
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(across(ends_with("_ImprintStatus"), ~ ifelse(. == "Paternally Imprinted", "P", .)))

# If value is "Maternally Imprinted" in columns ending in "_Imprinted_Status", change to "M"
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(across(ends_with("_ImprintStatus"), ~ ifelse(. == "Maternally Imprinted", "M", .)))


# Save the cleaned data frame as a new .csv file
write_csv(imprinting_analysis_clean, "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis_clean.csv")
