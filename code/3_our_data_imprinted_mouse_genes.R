###########################################################################
# Author: Kylee Duczyminski
# About: Finding all imprinted genes in the mouse genome and calculate the 
# percentage of mice that have paternal or maternal imprinting for each gene
# Input: data/processed/HIP_RNA-AlSp_parentalASE.csv
# Output: data/processed/HIP_RNA-AlSp_parentalASE_ourdata_imprinted.csv
###########################################################################

# Set proper working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyverse)

# Load in the data
imprinting_analysis <- read.csv("data/processed/HIP_RNA-AlSp_parentalASE.csv")

# Delete all columns ending in "_AER"
imprinting_analysis_clean <- imprinting_analysis %>%
  select(-ends_with("_AER"))

# If all values in a row are either "Biallelic" or "N/A", delete that row
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  filter(!(apply(.[, -1], 1, function(x) all(x == "Biallelic" | x == "N/A"))))

# Count the number of mice that are either "Maternally Imprinted" or "Paternally Imprinted" under columns ending in "_ImprintStatus" and create a column for this value
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  rowwise() %>%
  mutate(Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) %in% c("Maternally Imprinted", "Paternally Imprinted"), na.rm = TRUE)) %>%
  ungroup()

# Count the number of mice that have "Maternally Imprinted" under the columns with the ending "_ImprintStatus" and create a column for this value
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  rowwise() %>%
  mutate(Maternally_Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) == "Maternally Imprinted", na.rm = TRUE)) %>%
  ungroup()

# Count the number of mice that have "Paternally Imprinted" under the columns with the ending "_ImprintStatus" and create a column for this value
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  rowwise() %>%
  mutate(Paternally_Imprinted_Mice_Count = sum(c_across(ends_with("_ImprintStatus")) == "Paternally Imprinted", na.rm = TRUE)) %>%
  ungroup()

# Create a column called "Overall_Expression_Status" that is "Paternal" if all imprinted mice are maternally imprinted, "Maternal" if all imprinted mice are paternally imprinted, "Not Imprinted" if no mice are imprinted, and "Inconsistent" if there is a mix of maternally and paternally imprinted mice
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(Overall_Expression_Status = case_when(
    Maternally_Imprinted_Mice_Count == Imprinted_Mice_Count & Imprinted_Mice_Count > 0 ~ "Paternal",
    Paternally_Imprinted_Mice_Count == Imprinted_Mice_Count & Imprinted_Mice_Count > 0 ~ "Maternal",
    Imprinted_Mice_Count == 0 ~ "Not Imprinted",
    TRUE ~ "Inconsistent"
  ))

# Calculate the percent of imprinted mice that are maternally imprinted and create a column for this value
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(Percent_Maternally_Imprinted = ifelse(Imprinted_Mice_Count > 0, round((Maternally_Imprinted_Mice_Count / Imprinted_Mice_Count) * 100, 1), NA))

# Calculate the percent of imprinted mice that are paternally imprinted and create a column for this value
imprinting_analysis_clean <- imprinting_analysis_clean %>%
  mutate(Percent_Paternally_Imprinted = ifelse(Imprinted_Mice_Count > 0, round((Paternally_Imprinted_Mice_Count / Imprinted_Mice_Count) * 100, 1), NA))

# Save the cleaned data frame as a new .csv file
write_csv(imprinting_analysis_clean, "data/processed/HIP_RNA-AlSp_parentalASE_ourdata_imprinted.csv")
