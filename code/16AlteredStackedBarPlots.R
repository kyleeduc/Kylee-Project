############################################################################
# Author: Kylee Duczyminski
# Stacked barplot: WT observed imprinting consistency vs published expressed allele
# Input:
#   data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
#
# Output:
#   figures/WT_Only_Stacked_Barplot_Imprinting_Classification.jpg
############################################################################

# Load library
library(tidyverse)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Load data
imprinting_data <- read_csv(
  "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------------
# Use only the 6 WT samples
# ------------------------------------------------------------------
wt_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3"
)

wt_aer_columns <- paste0(wt_samples, "_AER")

# ------------------------------------------------------------------
# Clean published expressed allele
# ------------------------------------------------------------------
imprinting_data <- imprinting_data %>%
  mutate(
    Expressed_Allele = str_trim(Expressed_Allele),
    Expressed_Allele = case_when(
      Expressed_Allele %in% c("Maternal", "maternal", "Maternally Expressed") ~ "Maternal",
      Expressed_Allele %in% c("Paternal", "paternal", "Paternally Expressed") ~ "Paternal",
      TRUE ~ Expressed_Allele
    )
  ) %>%
  filter(Expressed_Allele %in% c("Maternal", "Paternal"))

# ------------------------------------------------------------------
# Function: convert AER into observed expressed allele
# ------------------------------------------------------------------
classify_observed_expression <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 0.7 ~ "Maternal",
    x <= -0.7 ~ "Paternal",
    TRUE ~ "Not Imprinted"
  )
}

# ------------------------------------------------------------------
# Long format: one row per gene per WT mouse sample
# ------------------------------------------------------------------
long_data <- imprinting_data %>%
  select(Gene_Name, Expressed_Allele, all_of(wt_aer_columns)) %>%
  pivot_longer(
    cols = all_of(wt_aer_columns),
    names_to = "Sample",
    values_to = "AER"
  ) %>%
  mutate(
    Sample = str_remove(Sample, "_AER"),
    Observed_Expressed_Allele = classify_observed_expression(AER)
  )

# ------------------------------------------------------------------
# Classify each gene based on all 6 WT samples
# ------------------------------------------------------------------
gene_classification <- long_data %>%
  group_by(Gene_Name, Expressed_Allele) %>%
  summarize(
    n_total = n(),
    n_maternal = sum(Observed_Expressed_Allele == "Maternal", na.rm = TRUE),
    n_paternal = sum(Observed_Expressed_Allele == "Paternal", na.rm = TRUE),
    n_not_imprinted = sum(Observed_Expressed_Allele == "Not Imprinted", na.rm = TRUE),
    n_imprinted = n_maternal + n_paternal,
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      
      # All 6 WT samples are not imprinted
      n_not_imprinted == n_total ~ "Not Imprinted",
      
      # All 6 WT samples match the published maternal allele
      Expressed_Allele == "Maternal" &
        n_maternal == n_total ~ "Consistent",
      
      # All 6 WT samples match the published paternal allele
      Expressed_Allele == "Paternal" &
        n_paternal == n_total ~ "Consistent",
      
      # All 6 WT samples are opposite of the published maternal allele
      Expressed_Allele == "Maternal" &
        n_paternal == n_total ~ "Opposite",
      
      # All 6 WT samples are opposite of the published paternal allele
      Expressed_Allele == "Paternal" &
        n_maternal == n_total ~ "Opposite",
      
      # Mixed maternal/paternal expression OR mixed imprinted/not imprinted calls
      (n_maternal > 0 & n_paternal > 0) |
        (n_imprinted > 0 & n_not_imprinted > 0) ~ "Inconsistent",
      
      # Catch anything unexpected
      TRUE ~ "Inconsistent"
    )
  )

# ------------------------------------------------------------------
# Summarize for stacked barplot
# ------------------------------------------------------------------
barplot_data <- gene_classification %>%
  mutate(
    Classification = str_trim(Classification),
    Expressed_Allele = str_trim(Expressed_Allele)
  ) %>%
  count(Expressed_Allele, Classification) %>%
  group_by(Expressed_Allele) %>%
  mutate(
    Proportion = n / sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    Classification = factor(
      Classification,
      levels = c("Consistent", "Opposite", "Inconsistent", "Not Imprinted")
    ),
    Expressed_Allele = factor(
      Expressed_Allele,
      levels = c("Maternal", "Paternal")
    )
  )

# ------------------------------------------------------------------
# Plot stacked barplot
# ------------------------------------------------------------------
stacked_barplot <- ggplot(
  barplot_data,
  aes(x = Expressed_Allele, y = Proportion, fill = Classification)
) +
  geom_col(color = "black", linewidth = 0.3) +
  
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  
  scale_fill_manual(
    values = c(
      "Consistent" = "#1E4E79",
      "Opposite" = "#FFCB05",
      "Inconsistent" = "#2E8B57",
      "Not Imprinted" = "#B0B0B0"
    ),
    drop = FALSE
  ) +
  
  labs(
    title = "WT Observed Expression Status Compared to Published\nExpressed Allele for Imprinted Genes",
    x = "Published Expressed Allele",
    y = "Proportion of Genes",
    fill = "Observed Classification"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

stacked_barplot

ggsave(
  "figures/WT_Only_Stacked_Barplot_Imprinting_Classification.jpg",
  stacked_barplot,
  width = 9,
  height = 5
)

##############################################################################
# Fisher's Exact Test
##############################################################################

# This tests whether maternally and paternally published genes differ
# in their proportion of Consistent vs Opposite classifications.
# Inconsistent and Not Imprinted genes are excluded from this test.

fisher_data <- gene_classification %>%
  filter(Classification %in% c("Consistent", "Opposite")) %>%
  count(Expressed_Allele, Classification) %>%
  complete(
    Expressed_Allele = c("Maternal", "Paternal"),
    Classification = c("Consistent", "Opposite"),
    fill = list(n = 0)
  )

fisher_table <- fisher_data %>%
  pivot_wider(
    names_from = Classification,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("Expressed_Allele") %>%
  as.matrix()

fisher_table

fisher_result <- fisher.test(fisher_table)

fisher_result
