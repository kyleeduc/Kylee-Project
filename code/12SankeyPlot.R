###############################################################################
# Author: Kylee Duczyminski
# About: Classify genes as imprinted/non-imprinted in our RNA-seq data using
# read thresholds, compare to published imprinting annotations, and generate
# a Sankey plot showing overlap.
###############################################################################

library(tidyverse)
library(readr)
library(ggalluvial)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# ============================================================================
# LOAD DATA
# ============================================================================

# Replace these paths/file names if needed
our_data <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")
published_data <- read_csv("data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv")

# ============================================================================
# IMPORTANT: ADJUST THESE COLUMN NAMES TO MATCH YOUR FILE
# ============================================================================

# Column containing gene names
gene_col <- "Gene_Name"

# Column containing published annotation
# Example values might include:
# "Maternally Imprinted", "Paternally Imprinted", "Not Imprinted", etc.
published_pattern_col <- "Overall_Expression_Pattern"

# Sample-specific maternal read columns
maternal_cols <- c(
  "WT_M1_Maternal", "WT_M2_Maternal", "WT_M3_Maternal",
  "WT_F1_Maternal", "WT_F2_Maternal", "WT_F3_Maternal",
  "HT_F1_Maternal", "HT_F2_Maternal", "HT_F3_Maternal",
  "KO_M1_Maternal", "KO_M2_Maternal", "KO_M3_Maternal"
)

# Sample-specific paternal read columns
paternal_cols <- c(
  "WT_M1_Paternal", "WT_M2_Paternal", "WT_M3_Paternal",
  "WT_F1_Paternal", "WT_F2_Paternal", "WT_F3_Paternal",
  "HT_F1_Paternal", "HT_F2_Paternal", "HT_F3_Paternal",
  "KO_M1_Paternal", "KO_M2_Paternal", "KO_M3_Paternal"
)

# ============================================================================
# CHECK COLUMN LENGTHS
# ============================================================================

if(length(maternal_cols) != length(paternal_cols)) {
  stop("maternal_cols and paternal_cols must have the same length.")
}

# ============================================================================
# SUBSET TO NEEDED COLUMNS
# ============================================================================

analysis_df <- our_data %>%
  select(all_of(c(gene_col, published_pattern_col, maternal_cols, paternal_cols)))

# ============================================================================
# CONVERT TO LONG FORMAT
# ============================================================================

sample_names <- str_remove(maternal_cols, "_Maternal$")

maternal_long <- analysis_df %>%
  select(all_of(c(gene_col, maternal_cols))) %>%
  pivot_longer(
    cols = all_of(maternal_cols),
    names_to = "Sample",
    values_to = "Maternal_Reads"
  ) %>%
  mutate(Sample = str_remove(Sample, "_Maternal$"))

paternal_long <- analysis_df %>%
  select(all_of(c(gene_col, paternal_cols))) %>%
  pivot_longer(
    cols = all_of(paternal_cols),
    names_to = "Sample",
    values_to = "Paternal_Reads"
  ) %>%
  mutate(Sample = str_remove(Sample, "_Paternal$"))

long_df <- maternal_long %>%
  inner_join(paternal_long, by = c(gene_col, "Sample")) %>%
  left_join(
    analysis_df %>% select(all_of(c(gene_col, published_pattern_col))),
    by = gene_col
  )

# ============================================================================
# CLASSIFY SAMPLE-LEVEL IMPRESSION STATUS
# Rules:
# - If total reads < 5 -> removed from analysis at sample level
# - If maternal/(maternal+paternal) >= 0.85 OR paternal/(maternal+paternal) >= 0.85
#   -> imprinted
# - Otherwise non-imprinted
# ============================================================================

long_df <- long_df %>%
  mutate(
    Maternal_Reads = as.numeric(Maternal_Reads),
    Paternal_Reads = as.numeric(Paternal_Reads),
    Total_Reads = Maternal_Reads + Paternal_Reads,
    Maternal_Fraction = ifelse(Total_Reads > 0, Maternal_Reads / Total_Reads, NA),
    Paternal_Fraction = ifelse(Total_Reads > 0, Paternal_Reads / Total_Reads, NA),
    Sample_Imprint_Status = case_when(
      is.na(Total_Reads) ~ NA_character_,
      Total_Reads < 5 ~ NA_character_,
      Maternal_Fraction >= 0.85 ~ "Imprinted",
      Paternal_Fraction >= 0.85 ~ "Imprinted",
      TRUE ~ "Non-imprinted"
    )
  )

# ============================================================================
# GENE-LEVEL CLASSIFICATION IN OUR DATA
# Rule:
# - If even one out of 12 valid samples is imprinted -> Imprinted in our data
# - Else -> Non-imprinted in our data
# - If no valid samples remain after filtering -> Removed
# ============================================================================

gene_summary <- long_df %>%
  group_by(.data[[gene_col]]) %>%
  summarize(
    Published_Pattern = first(.data[[published_pattern_col]]),
    Valid_Samples = sum(!is.na(Sample_Imprint_Status)),
    Imprinted_Samples = sum(Sample_Imprint_Status == "Imprinted", na.rm = TRUE),
    Our_Data_Status = case_when(
      Valid_Samples == 0 ~ "Removed",
      Imprinted_Samples >= 1 ~ "Imprinted",
      TRUE ~ "Non-imprinted"
    ),
    .groups = "drop"
  )

# ============================================================================
# COLLAPSE PUBLISHED PATTERN INTO IMP/NON-IMP
# You can edit these rules depending on your annotation wording
# ============================================================================

gene_summary <- gene_summary %>%
  mutate(
    Published_Status = case_when(
      str_to_lower(Published_Pattern) %in% c("not imprinted", "non-imprinted", "biallelic") ~ "Non-imprinted",
      is.na(Published_Pattern) ~ "Unknown",
      TRUE ~ "Imprinted"
    )
  )

# ============================================================================
# REMOVE GENES WITH UNKNOWN PUBLISHED STATUS OR KEEP THEM IF YOU WANT
# ============================================================================

gene_summary_plot <- gene_summary %>%
  filter(Published_Status %in% c("Imprinted", "Non-imprinted"),
         Our_Data_Status %in% c("Imprinted", "Non-imprinted"))

# ============================================================================
# COUNTS FOR OVERLAP
# ============================================================================

overlap_counts <- gene_summary_plot %>%
  count(Published_Status, Our_Data_Status)

print(overlap_counts)

# Helpful summary numbers
summary_numbers <- gene_summary_plot %>%
  summarize(
    Published_Imprinted = sum(Published_Status == "Imprinted"),
    Published_NonImprinted = sum(Published_Status == "Non-imprinted"),
    Our_Imprinted = sum(Our_Data_Status == "Imprinted"),
    Our_NonImprinted = sum(Our_Data_Status == "Non-imprinted"),
    Overlap_Imprinted_Imprinted = sum(Published_Status == "Imprinted" & Our_Data_Status == "Imprinted"),
    Overlap_Imprinted_NonImprinted = sum(Published_Status == "Imprinted" & Our_Data_Status == "Non-imprinted"),
    Overlap_NonImprinted_Imprinted = sum(Published_Status == "Non-imprinted" & Our_Data_Status == "Imprinted"),
    Overlap_NonImprinted_NonImprinted = sum(Published_Status == "Non-imprinted" & Our_Data_Status == "Non-imprinted")
  )

print(summary_numbers)

# ============================================================================
# SANKEY / ALUVIAL PLOT
# ============================================================================

sankey_df <- overlap_counts %>%
  rename(axis1 = Published_Status, axis2 = Our_Data_Status, Freq = n)

Sankey_plot <- ggplot(sankey_df,
                      aes(axis1 = axis1, axis2 = axis2, y = Freq)) +
  geom_alluvium(aes(fill = axis1), width = 0.25, alpha = 0.8) +
  geom_stratum(width = 0.25, fill = "gray90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Published Data", "Our RNA-seq Data"), expand = c(.1, .1)) +
  labs(
    title = "Overlap Between Published and RNA-seq Imprinting Classifications",
    subtitle = "Gene classified as imprinted in our data if ≥1 valid sample shows ≥85% reads from one parent",
    x = NULL,
    y = "Number of Genes",
    fill = "Published Classification"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

Sankey_plot

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

write_csv(gene_summary, "data/processed/gene_imprinting_summary_for_sankey.csv")
write_csv(overlap_counts, "data/processed/imprinting_overlap_counts.csv")

ggsave("figures/imprinting_sankey_plot.pdf", plot = Sankey_plot, width = 8, height = 5)
ggsave("figures/imprinting_sankey_plot.png", plot = Sankey_plot, width = 8, height = 5, dpi = 300)

