###############################################################################
# Collapse allele-specific counts before DESeq2
# Purpose:
#   Convert Gene_129 and Gene_JF1 rows into one total Gene row.
###############################################################################

library(tidyverse)
library(DESeq2)
library(ashr)

setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# ------------------------------------------------------------------
# 1. Read raw counts
# ------------------------------------------------------------------

raw_counts <- read_csv(
  "data/raw/HIP_RNA-AlSp_rawcounts.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------------
# 2. Collapse _129 and _JF1 rows into total gene-level counts
# ------------------------------------------------------------------

raw_counts_total <- raw_counts %>%
  mutate(
    Gene_Name = str_remove(Gene_Name, "_129$|_JF1$")
  ) %>%
  group_by(Gene_Name) %>%
  summarise(
    across(where(is.numeric), sum),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# 3. Save collapsed count table
# ------------------------------------------------------------------

write_csv(
  raw_counts_total,
  "data/processed/HIP_RNA_total_gene_counts_collapsed.csv"
)

# ------------------------------------------------------------------
# 4. Prepare count matrix
# ------------------------------------------------------------------

count_matrix <- raw_counts_total %>%
  column_to_rownames("Gene_Name") %>%
  as.matrix()

count_matrix <- round(count_matrix)

# ------------------------------------------------------------------
# 5. Create sample metadata
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Rename sequencing sample IDs to correct biological sample names
# ------------------------------------------------------------------

sample_name_map <- c(
  "1_S245"  = "WT_F1",
  "2_S246"  = "HT_F1",
  "3_S247"  = "HT_F2",
  "4_S248"  = "KO_M3",
  "5_S249"  = "WT_M1",
  "6_S250"  = "WT_M2",
  "7_S251"  = "HT_F3",
  "8_S252"  = "WT_F2",
  "9_S253"  = "WT_F3",
  "10_S254" = "KO_M1",
  "11_S255" = "WT_M3",
  "12_S256" = "KO_M2"
)

colnames(count_matrix) <- sample_name_map[colnames(count_matrix)]

# Reorder columns to match sample metadata
ordered_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3",
  "HT_F1", "HT_F2", "HT_F3",
  "KO_M1", "KO_M2", "KO_M3"
)

count_matrix <- count_matrix[, ordered_samples]

# ------------------------------------------------------------------
# 6. Run DESeq2
# ------------------------------------------------------------------

sample_metadata <- tibble(
  Sample = ordered_samples,
  Genotype = c(
    "WT_M", "WT_M", "WT_M",
    "WT_F", "WT_F", "WT_F",
    "HT_F", "HT_F", "HT_F",
    "KO_M", "KO_M", "KO_M"
  )
) %>%
  column_to_rownames("Sample")

# Check that count matrix columns and metadata rows match exactly
all(colnames(count_matrix) == rownames(sample_metadata))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ Genotype
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds$Genotype <- relevel(dds$Genotype, ref = "WT_M")

dds <- DESeq(dds)

# ------------------------------------------------------------------
# 7. Extract contrasts
# ------------------------------------------------------------------

res_KO_M_vs_WT_M <- results(
  dds,
  contrast = c("Genotype", "KO_M", "WT_M")
)

res_HT_F_vs_WT_F <- results(
  dds,
  contrast = c("Genotype", "HT_F", "WT_F")
)

res_KO_M_vs_HT_F <- results(
  dds,
  contrast = c("Genotype", "KO_M", "HT_F")
)

# ------------------------------------------------------------------
# 8. Shrink log2 fold changes
# ------------------------------------------------------------------

res_KO_M_vs_WT_M_shrunk <- lfcShrink(
  dds,
  contrast = c("Genotype", "KO_M", "WT_M"),
  res = res_KO_M_vs_WT_M,
  type = "ashr"
)

res_HT_F_vs_WT_F_shrunk <- lfcShrink(
  dds,
  contrast = c("Genotype", "HT_F", "WT_F"),
  res = res_HT_F_vs_WT_F,
  type = "ashr"
)

res_KO_M_vs_HT_F_shrunk <- lfcShrink(
  dds,
  contrast = c("Genotype", "KO_M", "HT_F"),
  res = res_KO_M_vs_HT_F,
  type = "ashr"
)

# ------------------------------------------------------------------
# 9. Convert results to tables
# ------------------------------------------------------------------

res_to_table <- function(res_obj) {
  as.data.frame(res_obj) %>%
    rownames_to_column("Gene_Name") %>%
    as_tibble() %>%
    arrange(padj)
}

KO_M_vs_WT_M_total <- res_to_table(res_KO_M_vs_WT_M_shrunk)
HT_F_vs_WT_F_total <- res_to_table(res_HT_F_vs_WT_F_shrunk)
KO_M_vs_HT_F_total <- res_to_table(res_KO_M_vs_HT_F_shrunk)

# ------------------------------------------------------------------
# 10. Save DESeq2 total gene-level results
# ------------------------------------------------------------------

write_csv(
  KO_M_vs_WT_M_total,
  "data/processed/DESeq2_TOTAL_KO_M_vs_WT_M_shrunk.csv"
)

write_csv(
  HT_F_vs_WT_F_total,
  "data/processed/DESeq2_TOTAL_HT_F_vs_WT_F_shrunk.csv"
)

write_csv(
  KO_M_vs_HT_F_total,
  "data/processed/DESeq2_TOTAL_KO_M_vs_HT_F_shrunk.csv"
)

###############################################################################
# End of script
###############################################################################