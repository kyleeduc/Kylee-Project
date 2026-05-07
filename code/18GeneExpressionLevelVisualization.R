##############################################################################
# Author: Kylee Duczyminski
# Purpose:
#   Filter raw allele-specific RNA-seq counts to published imprinted genes,
#   combine maternal and paternal allele counts into total expression counts,
#   and create a compact dot/strip plot showing expression across sample groups.
#
# *THIS SCRIPT CONTAINS RANDOM JITTER AND MAY NOT BE INFORMATIVE TO US*
# 
# Input files:
#   data/raw/HIP_RNA-AlSp_rawcounts.csv
#   data/raw/imprinted_gene_list.csv
#
# Output files:
#   data/processed/HIP_RNA-AlSp_rawcounts_imprinted_total_counts.csv
#   data/processed/HIP_RNA-AlSp_rawcounts_imprinted_total_counts_long.csv
#   figures/Imprinted_Genes_Total_RawCounts_DotPlot.pdf
#
# Notes:
#   In the raw counts file, each gene appears as two rows:
#       Gene_129 = maternal allele counts
#       Gene_JF1 = paternal allele counts
#
#   This script combines those two rows:
#       total gene expression = 129 counts + JF1 counts
#
#   Therefore, the final filtered count table should have one row per gene,
#   with a maximum of 112 genes from the published imprinted gene list.
##############################################################################

# ------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------
# 2. Set working directory
# ------------------------------------------------------------------
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# ------------------------------------------------------------------
# 3. Read input files
# ------------------------------------------------------------------
raw_counts <- read_csv(
  "data/raw/HIP_RNA-AlSp_rawcounts.csv",
  show_col_types = FALSE
)

imprinted_genes <- read_csv(
  "data/raw/imprinted_gene_list.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------------
# 4. Rename raw sequencing columns to biological sample names
# ------------------------------------------------------------------
# These are the original column names from the raw counts file.
old_names <- c(
  "1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250",
  "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256"
)

# These are the correct biological sample names.
new_names <- c(
  "WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2",
  "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2"
)

# Replace raw sequencing column names with biological sample names.
colnames(raw_counts)[match(old_names, colnames(raw_counts))] <- new_names

# ------------------------------------------------------------------
# 5. Define sample group information
# ------------------------------------------------------------------
sample_key <- tibble(
  Sample = c(
    "WT_M1", "WT_M2", "WT_M3",
    "WT_F1", "WT_F2", "WT_F3",
    "HT_F1", "HT_F2", "HT_F3",
    "KO_M1", "KO_M2", "KO_M3"
  ),
  Group = c(
    "WT Male", "WT Male", "WT Male",
    "WT Female", "WT Female", "WT Female",
    "HT Female", "HT Female", "HT Female",
    "KO Male", "KO Male", "KO Male"
  )
)

sample_columns <- sample_key$Sample

# ------------------------------------------------------------------
# 6. Clean published imprinted gene list
# ------------------------------------------------------------------
imprinted_genes <- imprinted_genes %>%
  mutate(
    Gene = str_trim(Gene),
    Status = str_trim(Status),
    Expressed_Allele = str_trim(Expressed_Allele)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

# ------------------------------------------------------------------
# 7. Rename first column in raw counts
# ------------------------------------------------------------------
# The first column contains gene names like:
#   Gene_129
#   Gene_JF1
names(raw_counts)[1] <- "Gene_Allele"

# ------------------------------------------------------------------
# 8. Parse gene name and allele source
# ------------------------------------------------------------------
raw_counts_parsed <- raw_counts %>%
  mutate(
    Gene_Allele = str_trim(Gene_Allele),
    
    Allele_Source = case_when(
      str_detect(Gene_Allele, "_129$") ~ "Maternal",
      str_detect(Gene_Allele, "_JF1$") ~ "Paternal",
      TRUE ~ NA_character_
    ),
    
    Gene = Gene_Allele %>%
      str_remove("_129$") %>%
      str_remove("_JF1$")
  )

# ------------------------------------------------------------------
# 9. Check whether parsing worked
# ------------------------------------------------------------------
cat("\nAllele source counts after parsing:\n")

raw_counts_parsed %>%
  count(Allele_Source) %>%
  print()

cat("\nFirst 10 parsed gene rows:\n")

raw_counts_parsed %>%
  select(Gene_Allele, Gene, Allele_Source) %>%
  head(10) %>%
  print()

# ------------------------------------------------------------------
# 10. Check gene overlap between raw counts and published imprinted list
# ------------------------------------------------------------------
raw_gene_names <- unique(raw_counts_parsed$Gene)
published_gene_names <- unique(imprinted_genes$Gene)

overlapping_genes <- intersect(raw_gene_names, published_gene_names)
missing_from_raw <- setdiff(published_gene_names, raw_gene_names)

cat("\nNumber of unique genes in raw counts:", length(raw_gene_names), "\n")
cat("Number of genes in published imprinted list:", length(published_gene_names), "\n")
cat("Number of overlapping genes:", length(overlapping_genes), "\n")

cat("\nPublished imprinted genes missing from raw counts:\n")
print(missing_from_raw)

# ------------------------------------------------------------------
# 11. Keep only published imprinted genes
# ------------------------------------------------------------------
filtered_counts <- raw_counts_parsed %>%
  filter(Gene %in% published_gene_names) %>%
  left_join(imprinted_genes, by = "Gene")

cat("\nFiltered allele-specific counts summary:\n")

filtered_counts %>%
  summarize(
    rows = n(),
    unique_genes = n_distinct(Gene),
    unique_gene_allele_rows = n_distinct(Gene_Allele)
  ) %>%
  print()

# ------------------------------------------------------------------
# 12. Combine maternal and paternal allele counts
# ------------------------------------------------------------------
# This collapses:
#   Gene_129
#   Gene_JF1
#
# into one row:
#   Gene = total counts from both alleles
filtered_counts_total <- filtered_counts %>%
  group_by(Gene, Status, Expressed_Allele) %>%
  summarize(
    across(
      all_of(sample_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("\nTotal-count filtered data summary:\n")

filtered_counts_total %>%
  summarize(
    rows = n(),
    unique_genes = n_distinct(Gene)
  ) %>%
  print()

# ------------------------------------------------------------------
# 13. Save one-row-per-gene total-count CSV
# ------------------------------------------------------------------
write_csv(
  filtered_counts_total,
  "data/processed/HIP_RNA-AlSp_rawcounts_imprinted_total_counts.csv"
)

# ------------------------------------------------------------------
# 14. Convert total counts to long format for plotting
# ------------------------------------------------------------------
long_counts_total <- filtered_counts_total %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "RawCount"
  ) %>%
  left_join(sample_key, by = "Sample") %>%
  mutate(
    RawCount = as.numeric(RawCount),
    log2_count = log2(RawCount + 1),
    
    Group = factor(
      Group,
      levels = c("WT Male", "WT Female", "HT Female", "KO Male")
    ),
    
    Expressed_Allele = factor(
      Expressed_Allele,
      levels = c("Maternal", "Paternal")
    )
  ) %>%
  filter(!is.na(Group))

cat("\nLong-format total-count data summary:\n")

long_counts_total %>%
  summarize(
    rows = n(),
    unique_genes = n_distinct(Gene),
    unique_samples = n_distinct(Sample),
    unique_groups = n_distinct(Group)
  ) %>%
  print()

cat("\nRows per sample group:\n")

long_counts_total %>%
  count(Group) %>%
  print()

# ------------------------------------------------------------------
# 15. Save long-format total-count CSV
# ------------------------------------------------------------------
write_csv(
  long_counts_total,
  "data/processed/HIP_RNA-AlSp_rawcounts_imprinted_total_counts_long.csv"
)

# ------------------------------------------------------------------
# 16. Order genes
# ------------------------------------------------------------------
# This keeps the gene order from the published imprinted gene list.
gene_order <- imprinted_genes %>%
  filter(Gene %in% unique(long_counts_total$Gene)) %>%
  pull(Gene) %>%
  unique()

long_counts_total <- long_counts_total %>%
  mutate(
    Gene = factor(Gene, levels = rev(gene_order))
  )

# ------------------------------------------------------------------
# 17. Calculate group means
# ------------------------------------------------------------------
group_means_total <- long_counts_total %>%
  group_by(Gene, Group, Expressed_Allele) %>%
  summarize(
    mean_raw = mean(RawCount, na.rm = TRUE),
    mean_log2 = mean(log2_count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nGroup means summary:\n")

group_means_total %>%
  summarize(
    rows = n(),
    unique_genes = n_distinct(Gene)
  ) %>%
  print()

# ------------------------------------------------------------------
# 18. Create dot/strip plot
# ------------------------------------------------------------------
# Each small dot = one replicate
# Large outlined dot = group mean
# Dot size = log2(total raw count + 1)
dot_plot_total <- ggplot(
  long_counts_total,
  aes(x = Group, y = Gene)
) +
  geom_point(
    aes(size = log2_count, color = Group),
    position = position_jitter(width = 0.16, height = 0.14),
    alpha = 0.75
  ) +
  
  geom_point(
    data = group_means_total,
    aes(x = Group, y = Gene, size = mean_log2),
    shape = 21,
    stroke = 0.7,
    fill = "white",
    color = "black",
    alpha = 0.4
  ) +
  
  facet_grid(
    Expressed_Allele ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_size_continuous(
    name = "log2(total count + 1)",
    range = c(0.8, 5)
  ) +
  
  labs(
    title = "Total Raw RNA-seq Counts for Published Imprinted Genes",
    subtitle = "Maternal and paternal allele counts combined; each small dot is one replicate",
    x = "Sample Group",
    y = "Gene",
    color = "Group"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

dot_plot_total

# ------------------------------------------------------------------
# 19. Save plot
# ------------------------------------------------------------------
ggsave(
  "figures/Imprinted_Genes_Total_RawCounts_DotPlot.jpg",
  dot_plot_total,
  width = 12,
  height = 18
)

##############################################################################
# End of script
##############################################################################