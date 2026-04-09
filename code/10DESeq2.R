###############################################################################
# Author: Kylee Duczyminski
# About: Code to run DESeq2 on the different genotypes and sexes in the RNA-seq
# dataset.
# Input: HIP_RNA-AlSp_rawcounts.csv from the data/raw folder
#        HIP_RNA-Percent_JF1.csv from the data/raw folder
# Output: DESeq2 results tables for each comparison in the data/processed folder
###############################################################################

# Load libraries
library(dplyr)
library(readr)
library(stringr)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# ============================================================================
# LOAD IN DATA
# ============================================================================

# Read in raw count data (gene expression counts)
raw_counts <- read_csv("data/raw/HIP_RNA-AlSp_rawcounts.csv")

# Read in percent JF1 data (allele-specific expression info)
percent_jf1 <- read_csv("data/raw/HIP_RNA-Percent_JF1.csv")

# ============================================================================
# CLEAN RAW COUNTS DATA
# ============================================================================

# Create Allele (129 or JF1) and Base_Gene (gene without suffix)
raw_counts_clean <- raw_counts %>%
  mutate(
    Allele = case_when(
      str_detect(Gene_Name, "_129$") ~ "129",
      str_detect(Gene_Name, "_JF1$") ~ "JF1",
      TRUE ~ NA_character_
    ),
    Base_Gene = str_remove(Gene_Name, "_(129|JF1)$")
  ) %>%
  filter(!is.na(Allele)) %>%  # keep only rows that have allele info
  distinct()                   # remove duplicate rows

# Keep only genes that have BOTH 129 and JF1 alleles
raw_counts_clean <- raw_counts_clean %>%
  group_by(Base_Gene) %>%
  filter(all(c("129", "JF1") %in% Allele)) %>%
  ungroup()

# ==========================================================================
# CLEAN PERCENT JF1 DATA
# ==========================================================================

# Remove duplicate gene rows
percent_jf1_clean <- percent_jf1 %>%
  distinct(Gene_Name, .keep_all = TRUE)

# ===========================================================================
# MERGE DATASETS
# ===========================================================================

# Join raw counts with percent JF1 using Base_Gene
merged_data <- raw_counts_clean %>%
  inner_join(percent_jf1_clean, by = c("Base_Gene" = "Gene_Name"))

# Check that each Base_Gene appears exactly twice (129 and JF1)
merged_data %>%
  count(Base_Gene) %>%
  filter(n != 2)

# Keep first 13 columns (Gene_Name + 12 samples)
merged_data <- merged_data %>%
  select(1:13)

# ==========================================================================
# RENAME SAMPLE COLUMNS
# ==========================================================================

# Original column names from dataset
old_names <- c("1_S245.x", "2_S246.x", "3_S247.x", "4_S248.x", "5_S249.x", "6_S250.x", 
               "7_S251.x", "8_S252.x", "9_S253.x", "10_S254.x", "11_S255.x", "12_S256.x")

new_names <- c("WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2", 
               "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2")

# Replace old column names with new ones
colnames(merged_data)[match(old_names, colnames(merged_data))] <- new_names

# ===========================================================================
# REORDER COLUMNS
# ===========================================================================

# Organize the columns by genotype and number: WT_M, WT_F, HT_F, KO_M
ordered_columns <- c("WT_M1", "WT_M2", "WT_M3", "WT_F1", "WT_F2", "WT_F3", 
                     "HT_F1", "HT_F2", "HT_F3", "KO_M1", "KO_M2", "KO_M3")

# Keep Gene_Name + ordered samples
merged_data <- merged_data %>% 
  select(Gene_Name, all_of(ordered_columns))

# ===========================================================================
# CREATE METADATA
# ===========================================================================

# Create metadata table using information from the merged_data
metadata <- data.frame(Sample = ordered_columns) %>%
  mutate(
    Genotype = str_extract(Sample, "^[^0-9]+"),   # WT_M, WT_F, HT_F, KO_M
    Sex = str_extract(Sample, "[MF]"),           # M or F
    
    Litter = c("C", "B", "B",
               "B", "B", "C",
               "B", "C", "B",
               "C", "B", "C"),
    Batch = c("A", "A", "B",
              "A", "B", "B",
              "A", "A", "B",
              "B", "B", "A")
  )

# Use the below information to correctly assign values to each sample

# sample_map <- data.frame(
#   Sample = c("WT_M1", "WT_M2", "WT_M3",
#              "WT_F1", "WT_F2", "WT_F3",
#              "HT_F1", "HT_F2", "HT_F3",
#              "KO_M1", "KO_M2", "KO_M3"),
#   
#   Sample_ID = c("12122-IV-5", "12122-IV-6", "12122-IV-11",
#                 "12122-IV-1", "12122-IV-8", "12122-IV-9",
#                 "12122-IV-2", "12122-IV-3", "12122-IV-7",
#                 "12122-IV-10", "12122-IV-12", "12122-IV-4"),
#   
#   Litter = c("C", "B", "B",
#              "B", "B", "C",
#              "B", "C", "B",
#              "C", "B", "C"),
#   
#   Batch = c("A", "A", "B",
#             "A", "B", "B",
#             "A", "A", "B",
#             "B", "B", "A")
# )

# ===========================================================================
# GENES DYSREGULATED IN KO/HT RELATIVE TO WT - GENOTYPE GROUPS CREATED
# ===========================================================================

# Collapse HT and KO into one group: "mut"
metadata$genotype <- ifelse(
  grepl("^WT", metadata$Genotype),
  "WT",
  "mut"
)

# Convert to factor (required for DESeq2)
metadata$genotype <- factor(metadata$genotype)

# ===========================================================================
# PREPARE COUNT MATRIX
# ===========================================================================

# Remove Gene_Name column to create numeric matrix
count_matrix <- merged_data[, -1]

# Set gene names as row names
rownames(count_matrix) <- merged_data$Gene_Name

# Convert to integer matrix (required for DESeq2)
count_matrix <- round(as.matrix(count_matrix))
mode(count_matrix) <- "integer"

# Set sample names as row names in metadata
rownames(metadata) <- metadata$Sample

# Check alignment between counts and metadata
all(colnames(count_matrix) == metadata$Sample)

# ===========================================================================
# DESeq2: WT vs Mut (ALL)
# ===========================================================================

# Load DESeq2 library
library(DESeq2)

# Adjusted p-value cutoff
PADJ <- 0.05

# Create DESeq2 dataset using genotype as the variable of interest
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ genotype
)

# Set WT as the default level for comparison
dds$genotype <- relevel(dds$genotype, ref = "WT")

# Run differential expression analysis
dds <- DESeq(dds)

# ==========================================================================
# EXTRACT RESULTS
# ==========================================================================

# Compare mutant vs WT
res <- results(dds, contrast = c("genotype", "mut", "WT"))

# Convert results to data frame
res_df <- as.data.frame(res)

# Add gene names as a column
res_df$Gene_Name <- rownames(res_df)

# ==========================================================================
# SAVE RESULTS
# ==========================================================================

# Save results to CSV file in the processed data folder
write_csv(res_df, "data/processed/DESeq2_WT_vs_Mutant.csv")

# ==========================================================================
# CREATE PCA PLOT
# ==========================================================================

## Make PCA Plot after DESeq2
library(ggplot2)

# Define colors for the two genotype groups
PCApalette <- c("WT" = "steelblue2", "mut" = "mediumorchid2")

# Variance stabilizing transformation
# This normalizes the count data for visualization like PCA
vsd <- vst(dds, blind = TRUE)

# Create PCA data using genotype as the grouping variable
pca_data <- plotPCA(vsd, intgroup = c("genotype"), returnData = TRUE)

# Save percent variance explained by PC1 and PC2
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Add metadata columns back in for plotting options
pca_data <- data.frame(
  pca_data,
  Sample = metadata$Sample,
  Sex = metadata$Sex,
  Batch = metadata$Batch,
  Litter = metadata$Litter
)

# Create PCA plot
PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color = genotype, shape = Sex)) +
  geom_point(size = 4.5) +
  ggtitle("PCA Plot of WT vs Mutant Samples") +
  scale_color_manual(values = PCApalette) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()

PCA_plot

# Save the PCA plot
ggsave("figures/WT_vs_Mutant_PCAplot.pdf", plot = PCA_plot, width = 5, height = 4)

# =========================
# RUN DESeq2 WITH GROUP DESIGN
# =========================

library(DESeq2)

# Make sure Genotype is a factor
metadata$Genotype <- factor(metadata$Genotype,
                            levels = c("WT_M", "KO_M", "WT_F", "HT_F"))

# Create DESeq2 dataset using group as the variable of interest
dds_group <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Genotype
)

# Run DESeq2
dds_group <- DESeq(dds_group)

# Check the available coefficient names
resultsNames(dds_group)


# =========================================================================
# EXTRACT RESULTS
# =========================================================================

# WT vs KO in males
KO_M_vs_WT_M <- results(dds_group, contrast = c("Genotype", "KO_M", "WT_M"))

# WT vs Het in females
HT_F_vs_WT_F <- results(dds_group, contrast = c("Genotype", "HT_F", "WT_F"))


# =========================================================================
# CONVERT RESULTS TO DATA FRAMES
# =========================================================================

# KO_M vs WT_M
KO_M_vs_WT_M_df <- as.data.frame(KO_M_vs_WT_M)
KO_M_vs_WT_M_df$Gene_Name <- rownames(KO_M_vs_WT_M_df)

# HT_F vs WT_F
HT_F_vs_WT_F_df <- as.data.frame(HT_F_vs_WT_F)
HT_F_vs_WT_F_df$Gene_Name <- rownames(HT_F_vs_WT_F_df)


# =========================================================================
# OPTIONAL: SORT BY ADJUSTED P VALUE
# =========================================================================

KO_M_vs_WT_M_df <- KO_M_vs_WT_M_df[order(KO_M_vs_WT_M_df$padj), ]
HT_F_vs_WT_F_df <- HT_F_vs_WT_F_df[order(HT_F_vs_WT_F_df$padj), ]


# =========================================================================
# OPTIONAL: KEEP ONLY SIGNIFICANT GENES
# =========================================================================

PADJ <- 0.05

KO_M_vs_WT_M_sig <- subset(KO_M_vs_WT_M_df, padj < PADJ)
HT_F_vs_WT_F_sig <- subset(HT_F_vs_WT_F_df, padj < PADJ)


# =========================================================================
# SAVE RESULTS
# =========================================================================

write_csv(KO_M_vs_WT_M_df, "data/processed/DESeq2_KO_M_vs_WT_M_all.csv")
write_csv(HT_F_vs_WT_F_df, "data/processed/DESeq2_HT_F_vs_WT_F_all.csv")

write_csv(KO_M_vs_WT_M_sig, "data/processed/DESeq2_KO_M_vs_WT_M_sig.csv")
write_csv(HT_F_vs_WT_F_sig, "data/processed/DESeq2_HT_F_vs_WT_F_sig.csv")

# =========================================================================
# SAVE NORMALIZED EXPRESSION TABLE
# =========================================================================

# Get normalized counts across all samples
normalized_counts <- counts(dds, normalized = TRUE)

# Convert to data frame and add gene names as a column
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$Gene_Name <- rownames(normalized_counts_df)

# Move Gene_Name to the first column
normalized_counts_df <- normalized_counts_df %>%
  select(Gene_Name, everything())

# Save comprehensive normalized expression table
write_csv(normalized_counts_df, "data/processed/Normalized_counts_all_samples.csv")


# =========================================================================
# LFC SHRINKAGE: WT vs MUT
# =========================================================================

# Check coefficient names first
resultsNames(dds)

# Shrink log2 fold changes for mutant vs WT
# Use the coefficient name that matches your resultsNames(dds) output
# May need to install - BiocManager::install("apeglm")
res_mut_vs_WT_shrunk <- lfcShrink(
  dds,
  coef = "genotype_mut_vs_WT",
  type = "apeglm"
)

# Convert to data frame and add gene names
res_mut_vs_WT_shrunk_df <- as.data.frame(res_mut_vs_WT_shrunk)
res_mut_vs_WT_shrunk_df$Gene_Name <- rownames(res_mut_vs_WT_shrunk_df)

# Save shrunk results
write_csv(res_mut_vs_WT_shrunk_df, "data/processed/DESeq2_mut_vs_WT_shrunk.csv")


# =========================================================================
# LFC SHRINKAGE: SEX-SPECIFIC COMPARISONS
# =========================================================================

# Check coefficient names for the Genotype model
resultsNames(dds_group)

# Shrink log2 fold changes for KO_M vs WT_M
# May need to install - BiocManager::install("ashr")
res_KO_M_vs_WT_M_shrunk <- lfcShrink(
  dds_group,
  contrast = c("Genotype", "KO_M", "WT_M"),
  type = "ashr"
)


# Shrink log2 fold changes for HT_F vs WT_F
res_HT_F_vs_WT_F_shrunk <- lfcShrink(
  dds_group,
  contrast = c("Genotype", "HT_F", "WT_F"),
  type = "ashr"
)

# Shrink log2 fold changes for KO_M vs HT_F
res_KO_M_vs_HT_F_shrunk <- lfcShrink(
  dds_group,
  contrast = c("Genotype", "KO_M", "HT_F"),
  type = "ashr"
)

# Convert each to data frame and add gene names
res_KO_M_vs_WT_M_shrunk_df <- as.data.frame(res_KO_M_vs_WT_M_shrunk)
res_KO_M_vs_WT_M_shrunk_df$Gene_Name <- rownames(res_KO_M_vs_WT_M_shrunk_df)

res_HT_F_vs_WT_F_shrunk_df <- as.data.frame(res_HT_F_vs_WT_F_shrunk)
res_HT_F_vs_WT_F_shrunk_df$Gene_Name <- rownames(res_HT_F_vs_WT_F_shrunk_df)

res_KO_M_vs_HT_F_shrunk_df <- as.data.frame(res_KO_M_vs_HT_F_shrunk)
res_KO_M_vs_HT_F_shrunk_df$Gene_Name <- rownames(res_KO_M_vs_HT_F_shrunk_df)

# Save shrunk results tables
write_csv(res_KO_M_vs_WT_M_shrunk_df, "data/processed/DESeq2_KO_M_vs_WT_M_shrunk.csv")
write_csv(res_HT_F_vs_WT_F_shrunk_df, "data/processed/DESeq2_HT_F_vs_WT_F_shrunk.csv")
write_csv(res_KO_M_vs_HT_F_shrunk_df, "data/processed/DESeq2_KO_M_vs_HT_F_shrunk.csv")


# =========================================================================
# SAVE UNSHRUNK RESULTS TABLES TOO
# =========================================================================

# WT vs mut
res_mut_vs_WT <- results(dds, contrast = c("genotype", "mut", "WT"))
res_mut_vs_WT_df <- as.data.frame(res_mut_vs_WT)
res_mut_vs_WT_df$Gene_Name <- rownames(res_mut_vs_WT_df)
write_csv(res_mut_vs_WT_df, "data/processed/DESeq2_mut_vs_WT_unshrunk.csv")

# KO_M vs WT_M
res_KO_M_vs_WT_M <- results(dds_group, contrast = c("Genotype", "KO_M", "WT_M"))
res_KO_M_vs_WT_M_df <- as.data.frame(res_KO_M_vs_WT_M)
res_KO_M_vs_WT_M_df$Gene_Name <- rownames(res_KO_M_vs_WT_M_df)
write_csv(res_KO_M_vs_WT_M_df, "data/processed/DESeq2_KO_M_vs_WT_M_unshrunk.csv")

# HT_F vs WT_F
res_HT_F_vs_WT_F <- results(dds_group, contrast = c("Genotype", "HT_F", "WT_F"))
res_HT_F_vs_WT_F_df <- as.data.frame(res_HT_F_vs_WT_F)
res_HT_F_vs_WT_F_df$Gene_Name <- rownames(res_HT_F_vs_WT_F_df)
write_csv(res_HT_F_vs_WT_F_df, "data/processed/DESeq2_HT_F_vs_WT_F_unshrunk.csv")

# KO_M vs HT_F
res_KO_M_vs_HT_F <- results(dds_group, contrast = c("Genotype", "KO_M", "HT_F"))
res_KO_M_vs_HT_F_df <- as.data.frame(res_KO_M_vs_HT_F)
res_KO_M_vs_HT_F_df$Gene_Name <- rownames(res_KO_M_vs_HT_F_df)
write_csv(res_KO_M_vs_HT_F_df, "data/processed/DESeq2_KO_M_vs_HT_F_unshrunk.csv")

# =========================================================================
# MA PLOTS
# =========================================================================

# MA plot for mut vs WT
pdf("figures/MA_mut_vs_WT.pdf", width = 6, height = 5)
plotMA(res_mut_vs_WT_shrunk, ylim = c(-5, 5), main = "MA Plot: Mut vs WT")
dev.off()

# MA plot for KO_M vs WT_M
pdf("figures/MA_KO_M_vs_WT_M.pdf", width = 6, height = 5)
plotMA(res_KO_M_vs_WT_M_shrunk, ylim = c(-5, 5), main = "MA Plot: KO_M vs WT_M")
dev.off()

# MA plot for HT_F vs WT_F
pdf("figures/MA_HT_F_vs_WT_F.pdf", width = 6, height = 5)
plotMA(res_HT_F_vs_WT_F_shrunk, ylim = c(-5, 5), main = "MA Plot: HT_F vs WT_F")
dev.off()

# MA plot for KO_M vs HT_F
pdf("figures/MA_KO_M_vs_HT_F.pdf", width = 6, height = 5)
plotMA(res_KO_M_vs_HT_F_shrunk, ylim = c(-5, 5), main = "MA Plot: KO_M vs HT_F")
dev.off()


# =========================================================================
# BEGIN HEATMAP OF SIGNIFICANT GENES
# =========================================================================

library(pheatmap)
library(tibble)

# Select significant genes from both KO_M_vs_WT_M and HT_F_vs_WT_F comparisons
sig_genes <- unique(c(
  KO_M_vs_WT_M_sig$Gene_Name,
  HT_F_vs_WT_F_sig$Gene_Name
))

# Subset normalized counts for significant genes
heatmap_data <- normalized_counts_df %>%
  filter(Gene_Name %in% sig_genes)

rownames(heatmap_data) <- heatmap_data$Gene_Name
heatmap_data$Gene_Name <- NULL

# Scale the data - Convert to Z-Scores
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Create annotation for samples
annotation_col <- data.frame(
  Genotype = metadata$Genotype
)

rownames(annotation_col) <- colnames(heatmap_data_scaled)

# Create and save heatmap
pdf("figures/heatmap_sig_genes.pdf", width = 8, height = 10)
pheatmap(
  heatmap_data_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none"
)

dev.off()

# =========================================================================
# CREATE LOG10 BASE MEAN MA PLOTS
# =========================================================================

library(ggplot2)

# MA plot for mut vs WT with log10 base mean on x-axis
ma_df <- as.data.frame(res_mut_vs_WT_shrunk)
ma_df$Gene_Name <- rownames(ma_df)
ma_df <- ma_df %>%
  mutate(log10_baseMean = log10(baseMean + 1))  # add 1 to avoid log(0)

ggplot(ma_df, aes(x = log10_baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  xlab("Log10(Base Mean + 1)") +
  ylab("Log2 Fold Change") +
  ggtitle("MA Plot: Mut vs WT with Log10 Base Mean")