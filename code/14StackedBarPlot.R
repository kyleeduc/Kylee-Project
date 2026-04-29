##############################################################################
# Stacked barplot: observed imprinting consistency vs published expressed allele
##############################################################################

library(tidyverse)

setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

imprinting_data <- read_csv(
  "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv",
  show_col_types = FALSE
)

ordered_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3",
  "HT_F1", "HT_F2", "HT_F3",
  "KO_M1", "KO_M2", "KO_M3"
)

ordered_aer_columns <- paste0(ordered_samples, "_AER")

# Clean published expressed allele
imprinting_data <- imprinting_data %>%
  mutate(
    Expressed_Allele = str_trim(Expressed_Allele),
    Expressed_Allele = case_when(
      Expressed_Allele %in% c("Maternal", "maternal", "Maternally Expressed") ~ "Maternal",
      Expressed_Allele %in% c("Paternal", "paternal", "Paternally Expressed") ~ "Paternal",
      TRUE ~ Expressed_Allele
    )
  )

# Function: convert AER into observed expressed allele
classify_observed_expression <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 0.7 ~ "Maternal",
    x <= -0.7 ~ "Paternal",
    TRUE ~ "Not Imprinted"
  )
}

# Long format: one row per gene per mouse sample
long_data <- imprinting_data %>%
  select(Gene_Name, Expressed_Allele, all_of(ordered_aer_columns)) %>%
  pivot_longer(
    cols = all_of(ordered_aer_columns),
    names_to = "Sample",
    values_to = "AER"
  ) %>%
  mutate(
    Sample = str_remove(Sample, "_AER"),
    Observed_Expressed_Allele = classify_observed_expression(AER)
  )

# Classify each gene as Consistent, Opposite, Inconsistent, or Not Imprinted
gene_classification <- long_data %>%
  group_by(Gene_Name, Expressed_Allele) %>%
  summarize(
    n_maternal = sum(Observed_Expressed_Allele == "Maternal", na.rm = TRUE),
    n_paternal = sum(Observed_Expressed_Allele == "Paternal", na.rm = TRUE),
    n_not_imprinted = sum(Observed_Expressed_Allele == "Not Imprinted", na.rm = TRUE),
    n_imprinted = n_maternal + n_paternal,
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      n_imprinted == 0 ~ "Not Imprinted",
      
      Expressed_Allele == "Maternal" & n_maternal > 0 & n_paternal == 0 ~ "Consistent",
      Expressed_Allele == "Paternal" & n_paternal > 0 & n_maternal == 0 ~ "Consistent",
      
      Expressed_Allele == "Maternal" & n_maternal == 0 & n_paternal > 0 ~ "Opposite",
      Expressed_Allele == "Paternal" & n_paternal == 0 & n_maternal > 0 ~ "Opposite",
      
      n_maternal > 0 & n_paternal > 0 ~ "Inconsistent",
      
      TRUE ~ "Inconsistent"
    )
  )

# Summarize for stacked barplot
barplot_data <- gene_classification %>%
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

# Plot stacked barplot
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
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c(
      "Consistent" = "#4daf4a",
      "Opposite" = "#e41a1c",
      "Inconsistent" = "#ff7f00",
      "Not Imprinted" = "gray70"
    )
  ) +
  labs(
    title = "Observed Imprinting Status Compared to Published Expressed Allele",
    x = "Published Expressed Allele",
    y = "Proportion of Genes",
    fill = "Observed Classification"
  ) +
  theme_classic(base_size = 14)

stacked_barplot

ggsave(
  "figures/Stacked_Barplot_Imprinting_Classification.pdf",
  stacked_barplot,
  width = 9,
  height = 5
)

##############################################################################
# Fisher's Exact Test
##############################################################################

# Collapse into matched vs mismatched/opposite only
fisher_data <- gene_classification %>%
  filter(Classification %in% c("Consistent", "Opposite")) %>%
  count(Expressed_Allele, Classification)

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
