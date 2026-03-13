##############################################################################
# Author: Kylee Duczyminski
# About: This script removes the rows of genes that are not imprinted
# Input: data/processed/HIP_RNA-mouse_imprinted_genes.csv
# Output: data/processed/mouse_imprinted_genes_list.csv
##############################################################################

# Load libraries
library(tidyverse)
library(readr)
library(dplyr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Read in the data
imprinted_genes <- read_csv("data/processed/HIP_RNA-mouse_imprinted_genes.csv")

# Remove the rows of genes that are not imprinted represented by the letter "N" in the "Imprinting_Category" column
imprinted_genes_list <- imprinted_genes %>%
  filter(Imprinting_Category != "N")

# Return the number of rows that remain in the file
num_rows <- nrow(imprinted_genes_list)
print(paste("Number of rows that remain in the file:", num_rows))

# Write the new file to the processed data folder
write_csv(imprinted_genes_list, "data/processed/mouse_imprinted_genes_list.csv")


