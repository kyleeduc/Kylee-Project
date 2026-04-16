#############################################################################
# Author: Kylee Duczyminski
# About: Which genes are dysregulated in Kdm5c -/Y and Kdm5c +/- mice compared 
# to Kdm5c +/Y and Kdm5c +/+ (General Dysregulation)
# Input: HIP_RNA-Counts.xlsx from the data/raw folder
#        HIP_RNA-Percent_JF1.csv from the data/raw folder
# Output: DESeq2 results for allele-specific expression in the data/processed folder
#############################################################################

# =========================================================================
# LOAD LIBRARIES
# =========================================================================

library(readr)
library(dplyr)
library(tidyr)

# =========================================================================
# SET WORKING DIRECTORY
# =========================================================================

setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# =========================================================================
# LOAD DATA
# =========================================================================

raw_counts <- read_csv("data/raw/HIP_RNA-AlSp_rawcounts.csv")
percent_jf1 <- read_csv("data/raw/HIP_RNA-Percent_JF1.csv")
snp_counts <- read_csv("data/raw/HIP_RNA-NumSNPs.csv")

# =========================================================================
# RENAME SAMPLE COLUMNS
# =========================================================================

old_names <- c("1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250",
               "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256")

new_names <- c("WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2",
               "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2")

colnames(percent_jf1)[match(old_names, colnames(percent_jf1))] <- new_names
colnames(snp_counts)[match(old_names, colnames(snp_counts))] <- new_names

# Put samples in the order you want
ordered_samples <- c("WT_M1", "WT_M2", "WT_M3",
                     "WT_F1", "WT_F2", "WT_F3",
                     "HT_F1", "HT_F2", "HT_F3",
                     "KO_M1", "KO_M2", "KO_M3")


