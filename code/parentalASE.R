##################################################################
# Author: Kylee Duczyminski
# About: Calculating parental allelic expression ratios in terms of paternal
# relative to maternal
# Input: HIP_RNA-AlSp_rawcounts.csv
# Output: HIP_RNA-AlSp_parentalASE.csv
##################################################################
library(dplyr)
library(tidyr)

# Set working directory to project folder
setwd("~/Documents/Iwase Lab/Kylee's Project")

# Load in raw counts .csv file
raw_counts <- read.csv("data/raw/HIP_RNA-AlSp_rawcounts.csv", header = TRUE)

# Maternal allelic expression has _129 at the end of Gene_Name
# Paternal allelic expression has _JF1 at the end of Gene_Name
# Create a new column with parental allelic expression ratios by calculating maternal minus paternal over maternal plus paternal

# IMPORTANT: Remove duplicate rows
raw_counts <- raw_counts %>%
  distinct()

# Step 1: Extract the base gene name and strain type
raw_counts <- raw_counts %>%
  mutate(
    base_gene_name = sub("_(129|JF1)$", "", Gene_Name),
    strain = ifelse(grepl("_129$", Gene_Name), "maternal", "paternal")
  )

# Step 2: Pivot to get maternal and paternal counts in separate columns
expression_wide <- raw_counts %>%
  select(-Gene_Name) %>%
  pivot_longer(
    cols = -c(base_gene_name, strain),
    names_to = "sample",
    values_to = "count"
  ) %>%
  pivot_wider(
    names_from = strain,
    values_from = count
  )

# Step 3: Calculate allelic expression ratio
expression_wide <- expression_wide %>%
  mutate(
    allelic_expression_ratio = (maternal - paternal) / (maternal + paternal),
    allelic_expression_ratio = ifelse(
      (maternal + paternal) == 0,
      NA,
      allelic_expression_ratio
    ),
    # Convert ratio to percentage (0 to 1 becomes 0 to 100)
    allelic_expression_ratio_percent = abs(allelic_expression_ratio) * 100,
    # Determine imprinting status (>85% absolute maternal or paternal bias)
    imprinting_status = ifelse(
      is.na(allelic_expression_ratio_percent),
      "Not Determined",
      ifelse(allelic_expression_ratio_percent > 85, "Imprinted", "Not Imprinted")
    )
  )

# Step 4: Pivot back to have genes as rows and samples as columns
results_final <- expression_wide %>%
  select(base_gene_name, sample, allelic_expression_ratio, allelic_expression_ratio_percent, imprinting_status) %>%
  pivot_wider(
    names_from = sample,
    values_from = c(allelic_expression_ratio, allelic_expression_ratio_percent, imprinting_status),
    names_glue = "{sample}_{.value}"
  )

# Step 5: View results
head(results_final)

# Step 6: Save to CSV
full_path <- "~/Documents/Iwase Lab/Kylee's Project/data/processed/HIP_RNA-AlSp_parentalASE.csv"
write.csv(results_final, full_path, row.names = FALSE)
cat("File saved to:", full_path, "\n")

# Verify the file was created
file.exists(full_path)

