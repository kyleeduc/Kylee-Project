##############################################################################
# Author: Kylee Duczyminski
# Stacked barplots: observed imprinting consistency vs published expressed allele
#
# Comparisons:
#   1. WT_M vs KO_M
#   2. WT_F vs HT_F
#   3. WT_M vs WT_F
#
# Classification is calculated separately within each genotype group.
#
# Categories:
#   Consistent     = all samples in that genotype match published allele
#   Opposite       = all samples in that genotype show opposite allele
#   Not Imprinted  = all samples in that genotype are not imprinted
#   Inconsistent   = mixed maternal/paternal OR mixed imprinted/not imprinted
##############################################################################

library(tidyverse)

setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

imprinting_data <- read_csv(
  "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv",
  show_col_types = FALSE
)

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
# Function: convert AER values into observed expressed allele calls
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
# Main function for each comparison
# ------------------------------------------------------------------
run_stacked_barplot_comparison <- function(sample_groups, group_name) {
  
  # sample_groups should be a named list, for example:
  # list(
  #   WT_M = c("WT_M1", "WT_M2", "WT_M3"),
  #   KO_M = c("KO_M1", "KO_M2", "KO_M3")
  # )
  
  selected_samples <- unlist(sample_groups)
  selected_aer_columns <- paste0(selected_samples, "_AER")
  
  sample_key <- tibble(
    Sample = selected_samples,
    Genotype = rep(names(sample_groups), lengths(sample_groups))
  )
  
  # --------------------------------------------------------------
  # Long format
  # --------------------------------------------------------------
  long_data <- imprinting_data %>%
    select(Gene_Name, Expressed_Allele, all_of(selected_aer_columns)) %>%
    pivot_longer(
      cols = all_of(selected_aer_columns),
      names_to = "Sample",
      values_to = "AER"
    ) %>%
    mutate(
      Sample = str_remove(Sample, "_AER"),
      Observed_Expressed_Allele = classify_observed_expression(AER)
    ) %>%
    left_join(sample_key, by = "Sample") %>%
    mutate(
      Genotype = factor(Genotype, levels = names(sample_groups))
    )
  
  # --------------------------------------------------------------
  # Classify each gene within each genotype group
  # --------------------------------------------------------------
  gene_classification <- long_data %>%
    group_by(Gene_Name, Expressed_Allele, Genotype) %>%
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
        
        # All samples in this genotype are not imprinted
        n_not_imprinted == n_total ~ "Not Imprinted",
        
        # All samples in this genotype match the published maternal allele
        Expressed_Allele == "Maternal" &
          n_maternal == n_total ~ "Consistent",
        
        # All samples in this genotype match the published paternal allele
        Expressed_Allele == "Paternal" &
          n_paternal == n_total ~ "Consistent",
        
        # All samples in this genotype are opposite of published maternal allele
        Expressed_Allele == "Maternal" &
          n_paternal == n_total ~ "Opposite",
        
        # All samples in this genotype are opposite of published paternal allele
        Expressed_Allele == "Paternal" &
          n_maternal == n_total ~ "Opposite",
        
        # Mixed maternal and paternal expression
        n_maternal > 0 & n_paternal > 0 ~ "Inconsistent",
        
        # Mixed imprinted and not imprinted calls
        n_imprinted > 0 & n_not_imprinted > 0 ~ "Inconsistent",
        
        # Catch anything unexpected
        TRUE ~ "Inconsistent"
      )
    )
  
  # ------------------------------------------------------------------
  # Read out genes classified as Opposite
  # ------------------------------------------------------------------
  opposite_genes <- gene_classification %>%
    filter(Classification == "Opposite") %>%
    select(
      Gene_Name,
      Genotype,
      Expressed_Allele,
      n_maternal,
      n_paternal,
      n_not_imprinted,
      n_imprinted,
      Classification
    ) %>%
    arrange(Genotype, Expressed_Allele, Gene_Name)
  
  cat("\n====================================================\n")
  cat(group_name, "Opposite Genes:\n")
  print(opposite_genes)
  cat("====================================================\n")
  
  write_csv(
    opposite_genes,
    paste0(
      "figures/",
      str_replace_all(group_name, "[^A-Za-z0-9]", "_"),
      "_Opposite_Genes.csv"
    )
  )
  
  # --------------------------------------------------------------
  # Summarize for stacked barplot
  # --------------------------------------------------------------
  barplot_data <- gene_classification %>%
    mutate(
      Classification = str_trim(Classification),
      Expressed_Allele = str_trim(Expressed_Allele)
    ) %>%
    count(Genotype, Expressed_Allele, Classification) %>%
    group_by(Genotype, Expressed_Allele) %>%
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
  
  # --------------------------------------------------------------
  # Make stacked barplot
  # --------------------------------------------------------------
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
    
    facet_wrap(~Genotype) +
    
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
      title = paste(
        group_name,
        "Observed Expression Status Compared to Published\nExpressed Allele for Imprinted Genes"
      ),
      x = "Published Expressed Allele",
      y = "Proportion of Genes",
      fill = "Observed Classification"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(stacked_barplot)
  
  # --------------------------------------------------------------
  # Save plot
  # --------------------------------------------------------------
  ggsave(
    filename = paste0(
      "figures/",
      str_replace_all(group_name, "[^A-Za-z0-9]", "_"),
      "_Stacked_Barplot.pdf"
    ),
    plot = stacked_barplot,
    width = 11,
    height = 5
  )
  
  # --------------------------------------------------------------
  # Fisher's Exact Test
  # --------------------------------------------------------------
  # Tests whether Consistent vs Opposite differs between genotype groups.
  # Inconsistent and Not Imprinted are excluded.
  
  fisher_data <- gene_classification %>%
    filter(Classification %in% c("Consistent", "Opposite")) %>%
    count(Genotype, Classification) %>%
    complete(
      Genotype,
      Classification = c("Consistent", "Opposite"),
      fill = list(n = 0)
    )
  
  fisher_table <- fisher_data %>%
    pivot_wider(
      names_from = Classification,
      values_from = n,
      values_fill = 0
    ) %>%
    column_to_rownames("Genotype") %>%
    as.matrix()
  
  fisher_result <- fisher.test(fisher_table)
  
  cat("\n====================================================\n")
  cat(group_name, "Fisher Table:\n")
  print(fisher_table)
  
  cat("\n", group_name, "Fisher's Exact Test:\n")
  print(fisher_result)
  cat("====================================================\n")
  
  return(
    list(
      long_data = long_data,
      gene_classification = gene_classification,
      barplot_data = barplot_data,
      plot = stacked_barplot,
      fisher_table = fisher_table,
      fisher_result = fisher_result
    )
  )
}

# ------------------------------------------------------------------
# 1. WT_M vs KO_M
# Classification based separately on WT_M1-3 and KO_M1-3
# ------------------------------------------------------------------
male_results <- run_stacked_barplot_comparison(
  sample_groups = list(
    WT_M = c("WT_M1", "WT_M2", "WT_M3"),
    KO_M = c("KO_M1", "KO_M2", "KO_M3")
  ),
  group_name = "WT_M_vs_KO_M"
)

# ------------------------------------------------------------------
# 2. WT_F vs HT_F
# Classification based separately on WT_F1-3 and HT_F1-3
# ------------------------------------------------------------------
female_results <- run_stacked_barplot_comparison(
  sample_groups = list(
    WT_F = c("WT_F1", "WT_F2", "WT_F3"),
    HT_F = c("HT_F1", "HT_F2", "HT_F3")
  ),
  group_name = "WT_F_vs_HT_F"
)

# ------------------------------------------------------------------
# 3. WT_M vs WT_F
# Classification based separately on WT_M1-3 and WT_F1-3
# ------------------------------------------------------------------
wt_sex_results <- run_stacked_barplot_comparison(
  sample_groups = list(
    WT_M = c("WT_M1", "WT_M2", "WT_M3"),
    WT_F = c("WT_F1", "WT_F2", "WT_F3")
  ),
  group_name = "WT_M_vs_WT_F"
)