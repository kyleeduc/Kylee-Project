###############################################################################
# Author: Kylee Duczyminski
# About: Compare how many of the 112 published imprinted genes have SNPs we can
# identify between 129 x JF1
# Input: data/raw/HIP_RNA-NumSNPs.csv
#        data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
# Output: data/processed/HIP_RNA-NumSNPsAnalysis.csv
###############################################################################

# Load libraries ----------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load in data
num_snps <- read_csv("data/raw/HIP_RNA-NumSNPs.csv")
gene_list <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# Merge datasets to compare the number of SNPs for imprinted vs non-imprinted genes
analysis_data <- gene_list %>%
  select(Gene_Name, Overall_Expression_Status, Data_Matching) %>%
  left_join(num_snps, by = "Gene_Name")

# Select sample columns
old_names <- c("1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250", 
               "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256")

new_names <- c("WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2", 
               "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2")

# Rename columns manually
colnames(analysis_data)[match(old_names, colnames(analysis_data))] <- new_names

# Organize the columns by genotype and number: WT_M, WT_F, HT_F, KO_M
ordered_columns <- c("WT_M1", "WT_M2", "WT_M3", "WT_F1", "WT_F2", "WT_F3", 
                     "HT_F1", "HT_F2", "HT_F3", "KO_M1", "KO_M2", "KO_M3")
analysis_data <- analysis_data %>% select(Gene_Name, Overall_Expression_Status, Data_Matching, all_of(ordered_columns))

# Select only the sample columns (exclude Gene_Name)
sample_columns <- analysis_data %>% select(all_of(ordered_columns))

# Convert all NA values to 0 in analysis_data
analysis_data[is.na(analysis_data)] <- 0

# Sum the number of SNPs for each gene across all samples
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(Total_SNPs = sum(c_across(all_of(ordered_columns))))

# Average the number of SNPs that were used across all 12 samples for each gene
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(Avg_SNPs = round(mean(c_across(all_of(ordered_columns)))))

# Sort the rows in ascending order from lowest Avg_SNPs to highest Avg_SNPs
analysis_data <- analysis_data %>%
  arrange(Total_SNPs)

# Generate one histogram that shows the average SNPs (y-axis), binned in groups of 10 SNPs (i.e., 20-29, 30-39) and number of genes (x-axis) for both conserved and nonconserved dataframes (use different colors for each, make the colors transparent so we can see the overlap). Save as a figure. [FIG1]
# Create bins for average SNPs
analysis_data <- analysis_data %>%
  mutate(SNP_Bin = cut(Avg_SNPs, breaks = seq(0, max(Avg_SNPs) + 10, by = 10), right = FALSE))

# Count the number of genes in each SNP bin for conserved and non-conserved genes. Conserved if Data_Matching == "Yes", non-conserved if Data_Matching == "No" or "Inconsistent"
bin_counts <- analysis_data %>%
  mutate(Conservation_Status = case_when(
    Data_Matching == "Yes" ~ "Conserved",
    Data_Matching %in% c("No", "Inconsistent") ~ "Non-Conserved",
    TRUE ~ NA_character_
  )) %>%
  group_by(SNP_Bin, Conservation_Status) %>%
  summarise(Gene_Count = n(), .groups = 'drop') %>%
  # Create all combinations of SNP_Bin and Conservation_Status
  complete(SNP_Bin, Conservation_Status, fill = list(Gene_Count = 0))


# Create a histogram
library(ggplot2)

ggplot(bin_counts, aes(x = SNP_Bin, y = Gene_Count, fill = Conservation_Status)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  labs(title = "Distribution of Average SNPs per Gene by Conservation Status",
       x = "Average Number of SNPs (binned)",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Conserved" = "blue", "Non-Conserved" = "red")) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2)
  )

# Save this histogram as a .pdf
ggsave("figures/HIP_RNA-NumSNPs_Histogram.pdf", width = 8, height = 6)

# Save the analysis_data dataframe as a .csv file in the processed data folder
write_csv(analysis_data, "data/processed/HIP_RNA-NumSNPsAnalysis.csv")

