##############################################################################
# Author: Kylee Duczyminski
# Purpose:
# Create separate stacked barplots and Fisher's Exact Tests for:
#   1. Male mice: WT_M vs KO_M
#   2. Female mice: WT_F vs HT_F
#   3. Wild-type sex comparison: WT_M vs WT_F
#
# The plot compares observed allele-specific expression patterns to the
# published expressed allele for each imprinted gene.
#
# Categories:
#   Consistent     = all samples match the published expressed allele
#   Opposite       = all observed expressed alleles are opposite of published data
#   Inconsistent   = some samples are maternal, some are paternal OR some samples are imprinted and some are not
#   Not Imprinted  = no samples pass the AER imprinting threshold
#
# Input:
#   data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
#
# Output:
#   figures/Male_WT_vs_KO_Stacked_Barplot.pdf
#   figures/Female_WT_vs_HT_Stacked_Barplot.pdf
#   figures/WT_M_vs_WT_F_Stacked_Barplot.pdf
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
# 3. Read in data
# ------------------------------------------------------------------
imprinting_data <- read_csv(
  "data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------------
# 4. Clean published expressed allele labels
# ------------------------------------------------------------------
# This makes sure labels are consistent as either "Maternal" or "Paternal".
imprinting_data <- imprinting_data %>%
  mutate(
    Expressed_Allele = str_trim(Expressed_Allele),
    Expressed_Allele = case_when(
      Expressed_Allele %in% c("Maternal", "maternal", "Maternally Expressed") ~ "Maternal",
      Expressed_Allele %in% c("Paternal", "paternal", "Paternally Expressed") ~ "Paternal",
      TRUE ~ Expressed_Allele
    )
  )

# ------------------------------------------------------------------
# 5. Define sample groups
# ------------------------------------------------------------------
# Male comparison: WT male vs KO male
male_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "KO_M1", "KO_M2", "KO_M3"
)

# Female comparison: WT female vs heterozygous female
female_samples <- c(
  "WT_F1", "WT_F2", "WT_F3",
  "HT_F1", "HT_F2", "HT_F3"
)

# Add "_AER" to match column names in the dataset
male_aer <- paste0(male_samples, "_AER")
female_aer <- paste0(female_samples, "_AER")

# WT comparisons: WT_M vs WT_F
wt_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3"
)

# Add "_AER" to match column names in the dataset
wt_aer <- paste0(wt_samples, "_AER")

# ------------------------------------------------------------------
# 6. Convert AER values into observed expressed allele calls
# ------------------------------------------------------------------
# AER >= 0.7  = maternally expressed
# AER <= -0.7 = paternally expressed
# Between -0.7 and 0.7 = not imprinted
classify_observed_expression <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 0.7 ~ "Maternal",
    x <= -0.7 ~ "Paternal",
    TRUE ~ "Not Imprinted"
  )
}

# ------------------------------------------------------------------
# 7. Main analysis function
# ------------------------------------------------------------------
# This function:
#   1. Converts data to long format
#   2. Classifies each gene within each genotype group
#   3. Makes a stacked barplot
#   4. Saves the plot as a PDF
#   5. Runs Fisher's Exact Test
run_group_analysis <- function(aer_columns, group_name) {
  
  # --------------------------------------------------------------
  # Convert wide AER data into long format
  # --------------------------------------------------------------
  long_data <- imprinting_data %>%
    select(Gene_Name, Expressed_Allele, all_of(aer_columns)) %>%
    pivot_longer(
      cols = all_of(aer_columns),
      names_to = "Sample",
      values_to = "AER"
    ) %>%
    mutate(
      Sample = str_remove(Sample, "_AER"),
      Observed_Expressed_Allele = classify_observed_expression(AER),
      
      # Assign genotype labels
      Genotype = case_when(
        str_detect(Sample, "WT_M") ~ "WT_M",
        str_detect(Sample, "WT_F") ~ "WT_F",
        str_detect(Sample, "KO_M") ~ "KO_M",
        str_detect(Sample, "HT_F") ~ "HT_F",
        TRUE ~ NA_character_
      )
    )
  
  if (group_name == "Male_WT_vs_KO") {
    long_data$Genotype <- factor(long_data$Genotype, levels = c("WT_M", "KO_M"))
  }
  
  if (group_name == "Female_WT_vs_HT") {
    long_data$Genotype <- factor(long_data$Genotype, levels = c("WT_F", "HT_F"))
  }
  
  if (group_name == "WT_M_vs_WT_F") {
    long_data$Genotype <- factor(long_data$Genotype, levels = c("WT_M", "WT_F"))
  }
  
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
        
        # No samples pass the imprinting threshold
        n_imprinted == 0 ~ "Not Imprinted",
        
        # All samples match the published maternal allele
        Expressed_Allele == "Maternal" &
          n_maternal == n_total ~ "Consistent",
        
        # All samples match the published paternal allele
        Expressed_Allele == "Paternal" &
          n_paternal == n_total ~ "Consistent",
        
        # All samples are opposite of published maternal allele
        Expressed_Allele == "Maternal" &
          n_paternal == n_total ~ "Opposite",
        
        # All samples are opposite of published paternal allele
        Expressed_Allele == "Paternal" &
          n_maternal == n_total ~ "Opposite",
        
        # Some samples are maternal and some are paternal
        n_maternal > 0 & n_paternal > 0 ~ "Inconsistent",
        
        # Some samples are imprinted and some are not imprinted
        n_imprinted > 0 & n_not_imprinted > 0 ~ "Inconsistent",
        
        # Catch anything unexpected
        TRUE ~ "Inconsistent"
      )
    )

  # --------------------------------------------------------------
  # Summarize data for stacked barplot
  # --------------------------------------------------------------
  barplot_data <- gene_classification %>%
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
  p <- ggplot(
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
      )
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
  
  print(p)
  
  # --------------------------------------------------------------
  # Save plot as PDF
  # --------------------------------------------------------------
  ggsave(
    filename = paste0(
      "figures/",
      str_replace_all(group_name, "[^A-Za-z0-9]", "_"),
      "_Stacked_Barplot.pdf"
    ),
    plot = p,
    width = 11,
    height = 8
  )
  
  # --------------------------------------------------------------
  # Fisher's Exact Test
  # --------------------------------------------------------------
  # This tests whether the proportion of Consistent vs Opposite genes
  # differs between genotype groups.
  #
  # Inconsistent and Not Imprinted genes are excluded from this test.
  
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
      plot = p,
      fisher_table = fisher_table,
      fisher_result = fisher_result
    )
  )
}


# ------------------------------------------------------------------
# 8. Run male comparison: WT_M vs KO_M
# ------------------------------------------------------------------
male_results <- run_group_analysis(
  male_aer,
  "Male_WT_vs_KO"
)

# ------------------------------------------------------------------
# 9. Run female comparison: WT_F vs HT_F
# ------------------------------------------------------------------
female_results <- run_group_analysis(
  female_aer,
  "Female_WT_vs_HT"
)

# ------------------------------------------------------------------
# 10. Run WT comparison: WT_M vs WT_F
# ------------------------------------------------------------------
wt_results <- run_group_analysis(
  wt_aer,
  "WT_M_vs_WT_F"
)
