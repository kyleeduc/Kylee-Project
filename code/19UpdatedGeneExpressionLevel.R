##############################################################################
# Author: Kylee Duczyminski
# Purpose:
#   Filter raw allele-specific RNA-seq counts to published imprinted genes,
#   calculate allelic expression ratio (AER), and create a dot plot showing
#   whether each sample is maternally expressed, paternally expressed, or not
#   imprinted.
#
# Input files:
#   data/raw/HIP_RNA-AlSp_rawcounts.csv
#   data/raw/imprinted_gene_list.csv
#
# Output files:
#   data/processed/HIP_RNA-AlSp_rawcounts_imprinted_AER_long.csv
#   figures/Imprinted_Genes_AER_DotPlot.jpg
#
# Notes:
#   In the raw counts file, each gene appears as two rows:
#       Gene_129 = maternal allele counts
#       Gene_JF1 = paternal allele counts
#
#   AER is calculated as:
#       AER = (Maternal_Count - Paternal_Count) /
#             (Maternal_Count + Paternal_Count)
#
#   Interpretation:
#       AER close to +1 = mostly maternal / 129 expression
#       AER close to -1 = mostly paternal / JF1 expression
#       AER close to  0 = balanced / not strongly imprinted
#
#   Classification thresholds:
#       AER >=  0.7 = Maternal
#       AER <= -0.7 = Paternal
#       Otherwise   = Not Imprinted
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
# Original column names from the raw counts file
old_names <- c(
  "1_S245", "2_S246", "3_S247", "4_S248", "5_S249", "6_S250",
  "7_S251", "8_S252", "9_S253", "10_S254", "11_S255", "12_S256"
)

# Correct biological sample names
new_names <- c(
  "WT_F1", "HT_F1", "HT_F2", "KO_M3", "WT_M1", "WT_M2",
  "HT_F3", "WT_F2", "WT_F3", "KO_M1", "WT_M3", "KO_M2"
)

# Replace raw sequencing column names with biological sample names
colnames(raw_counts)[match(old_names, colnames(raw_counts))] <- new_names

# ------------------------------------------------------------------
# 5. Define sample groups
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
    Expressed_Allele = str_trim(Expressed_Allele),
    Expressed_Allele = case_when(
      Expressed_Allele %in% c("Maternal", "maternal", "Maternally Expressed") ~ "Maternal",
      Expressed_Allele %in% c("Paternal", "paternal", "Paternally Expressed") ~ "Paternal",
      TRUE ~ Expressed_Allele
    )
  ) %>%
  distinct(Gene, .keep_all = TRUE)

# ------------------------------------------------------------------
# 7. Rename first column in raw counts
# ------------------------------------------------------------------
# First column contains values like:
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
# 12. Convert to long format
# ------------------------------------------------------------------
# This keeps allele source separate at first.
long_allele_counts <- filtered_counts %>%
  select(
    Gene,
    Gene_Allele,
    Allele_Source,
    Status,
    Expressed_Allele,
    all_of(sample_columns)
  ) %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "RawCount"
  ) %>%
  left_join(sample_key, by = "Sample") %>%
  mutate(
    RawCount = as.numeric(RawCount)
  ) %>%
  filter(!is.na(Group))

# ------------------------------------------------------------------
# 13. Convert maternal and paternal rows into columns
# ------------------------------------------------------------------
# This creates one row per gene per sample:
#   Maternal = 129 counts
#   Paternal = JF1 counts
aer_long <- long_allele_counts %>%
  select(
    Gene,
    Status,
    Expressed_Allele,
    Sample,
    Group,
    Allele_Source,
    RawCount
  ) %>%
  pivot_wider(
    names_from = Allele_Source,
    values_from = RawCount,
    values_fill = 0
  )

# ------------------------------------------------------------------
# 14. Calculate total counts and AER
# ------------------------------------------------------------------
aer_long <- aer_long %>%
  mutate(
    Total_Count = Maternal + Paternal,
    
    AER = case_when(
      Total_Count == 0 ~ NA_real_,
      TRUE ~ (Maternal - Paternal) / Total_Count
    ),
    
    Observed_Expression = case_when(
      Total_Count == 0 ~ "No Expression",
      AER >= 0.7 ~ "Maternal",
      AER <= -0.7 ~ "Paternal",
      TRUE ~ "Not Imprinted"
    ),
    
    log2_total_count = log2(Total_Count + 1),
    
    Group = factor(
      Group,
      levels = c("WT Male", "WT Female", "HT Female", "KO Male")
    ),
    
    Expressed_Allele = factor(
      Expressed_Allele,
      levels = c("Maternal", "Paternal")
    ),
    
    Observed_Expression = factor(
      Observed_Expression,
      levels = c("Maternal", "Paternal", "Not Imprinted", "No Expression")
    )
  )

# ------------------------------------------------------------------
# 15. Save AER long-format file
# ------------------------------------------------------------------
write_csv(
  aer_long,
  "data/processed/HIP_RNA-AlSp_rawcounts_imprinted_AER_long.csv"
)

# ------------------------------------------------------------------
# 16. Order genes by average total expression
# ------------------------------------------------------------------
# Highest-expressed genes will appear at the TOP of the plot.
# We use the average log2(total count + 1) across all samples.

gene_expression_order <- aer_long %>%
  group_by(Gene) %>%
  summarize(
    mean_log2_total = mean(log2_total_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_log2_total))

# In ggplot, the first factor level appears at the bottom,
# so we reverse the order so the highest-expression genes are at the top.
aer_long <- aer_long %>%
  mutate(
    Gene = factor(Gene, levels = rev(gene_expression_order$Gene))
  )

# ------------------------------------------------------------------
# 17. Keep only rows with valid AER for plotting
# ------------------------------------------------------------------
# No-expression rows cannot be plotted on an AER axis because AER is NA.
aer_plot_data <- aer_long %>%
  filter(!is.na(AER))

# ------------------------------------------------------------------
# 18. Calculate group means for optional summary dots
# ------------------------------------------------------------------
group_means_aer <- aer_plot_data %>%
  group_by(Gene, Group, Expressed_Allele) %>%
  summarize(
    mean_AER = mean(AER, na.rm = TRUE),
    mean_total_count = mean(Total_Count, na.rm = TRUE),
    mean_log2_total_count = mean(log2_total_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = levels(aer_plot_data$Gene))
  )

# ------------------------------------------------------------------
# 19. Create AER dot plot
# ------------------------------------------------------------------
# Interpretation:
#   Left side  = paternal / JF1 expression
#   Right side = maternal / 129 expression
#
# Dot color:
#   Navy  = Maternal
#   Maize = Paternal
#   Gray  = Not Imprinted
#   Light gray = No Expression (not shown on AER axis because AER is NA)
#
# Dot size:
#   Larger dots = higher total expression

aer_dot_plot <- ggplot(
  aer_plot_data,
  aes(x = AER, y = Gene)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = c(-0.7, 0.7),
    linetype = "dotted",
    color = "gray40",
    linewidth = 0.4
  ) +
  
  geom_point(
    aes(
      color = Observed_Expression,
      size = log2_total_count
    ),
    alpha = 0.75,
    position = position_jitter(width = 0.03, height = 0.14)
  ) +
  
  # Mean point for each gene/group
  geom_point(
    data = group_means_aer,
    aes(
      x = mean_AER,
      y = Gene,
      size = mean_log2_total_count
    ),
    shape = 21,
    stroke = 0.7,
    fill = "white",
    color = "black",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~Group, nrow = 1) +
  
  scale_x_continuous(
    breaks = c(-1, -0.7, 0, 0.7, 1),
    labels = c(
      "Paternal\n-1",
      "-0.7",
      "Balanced\n0",
      "0.7",
      "Maternal\n+1"
    )
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  
  scale_color_manual(
    values = c(
      "Maternal" = "#FFCB05",
      "Paternal" = "#1E4E79",
      "Not Imprinted" = "#B0B0B0",
      "No Expression" = "#E6E6E6"
    ),
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name = "log2(total count + 1)",
    range = c(1.2, 7)
  ) +
  
  labs(
    title = "Allele-Specific Expression Bias for Published Imprinted Genes",
    subtitle = "Genes ordered by average total expression (highest at top)\nLeft = paternal/JF1 expression; right = maternal/129 expression",
    x = "Allelic Expression Ratio",
    y = "Gene",
    color = "Observed Expression"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.spacing.x = grid::unit(1.4, "lines")
  )

aer_dot_plot

# ------------------------------------------------------------------
# 20. Save AER dot plot
# ------------------------------------------------------------------
ggsave(
  "figures/Imprinted_Genes_AER_DotPlot_Ordered.pdf",
  aer_dot_plot,
  width = 18,
  height = 20,
  dpi = 300
)

##############################################################################
# End of script
##############################################################################