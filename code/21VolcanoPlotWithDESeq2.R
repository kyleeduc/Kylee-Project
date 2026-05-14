###############################################################################
# Author: Kylee Duczyminski
# Purpose:
#   Create volcano plots and violin plots from saved DESeq2 results.
#
# Volcano plot rules:
#   - x-threshold: |log2FC| > 0.25
#   - y-threshold: padj < 0.05
#
# Point coloring:
#   - Significant imprinted maternal genes = deeppink3
#   - Significant imprinted paternal genes = dodgerblue3
#   - Significant non-imprinted DEGs = grey25
#   - Small fold change genes = grey25
#   - Not significant genes = lightgray
#
# Labels:
#   - Label only significant imprinted DEGs
#   - Label text = Gene_Name
#
# Inputs:
#   data/processed/DESeq2_KO_M_vs_WT_M_shrunk.csv
#   data/processed/DESeq2_HT_F_vs_WT_F_shrunk.csv
#   data/processed/DESeq2_KO_M_vs_HT_F_shrunk.csv
#   data/processed/Normalized_counts_all_samples.csv
#   data/raw/imprinted_gene_list.csv
#
# Outputs:
#   figures/Volcano_KO_M_vs_WT_M.pdf
#   figures/Volcano_HT_F_vs_WT_F.pdf
#   figures/Volcano_KO_M_vs_HT_F.pdf
#   figures/Violin_Imprinted_DEGs_KO_M_vs_WT_M.pdf
#   figures/Violin_Imprinted_DEGs_HT_F_vs_WT_F.pdf
#   figures/Violin_Imprinted_DEGs_KO_M_vs_HT_F.pdf
###############################################################################

library(tidyverse)
library(ggrepel)
library(ggplot2)

setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# ------------------------------------------------------------------
# 1. Set thresholds and colors
# ------------------------------------------------------------------

PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 0.25

MATERNAL_COLOR <- "deeppink3"
PATERNAL_COLOR <- "dodgerblue3"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 2. Read published imprinted gene list
# ------------------------------------------------------------------

imprinted_genes_raw <- read_csv(
  "data/raw/imprinted_gene_list.csv",
  show_col_types = FALSE
)

imprinted_genes <- imprinted_genes_raw %>%
  mutate(
    Gene = str_trim(Gene)
  )

if (!"Expressed_Allele" %in% names(imprinted_genes)) {
  imprinted_genes$Expressed_Allele <- NA_character_
}

if (!"Overall_Expression_Pattern" %in% names(imprinted_genes)) {
  imprinted_genes$Overall_Expression_Pattern <- NA_character_
}

imprinted_genes <- imprinted_genes %>%
  mutate(
    Expressed_Allele = str_trim(Expressed_Allele),
    Overall_Expression_Pattern = str_trim(Overall_Expression_Pattern),
    
    Expressed_Allele = case_when(
      Expressed_Allele %in% c(
        "Maternal", "maternal", "Maternally Expressed",
        "Maternally expressed", "Maternal expression"
      ) ~ "Maternal",
      
      Expressed_Allele %in% c(
        "Paternal", "paternal", "Paternally Expressed",
        "Paternally expressed", "Paternal expression"
      ) ~ "Paternal",
      
      TRUE ~ Expressed_Allele
    )
  ) %>%
  filter(
    Expressed_Allele %in% c("Maternal", "Paternal")
  ) %>%
  filter(
    is.na(Overall_Expression_Pattern) |
      Overall_Expression_Pattern != "Not Imprinted"
  ) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  select(Gene, Expressed_Allele, everything())

# ------------------------------------------------------------------
# 3. Read DESeq2 result tables
# ------------------------------------------------------------------

KO_M_vs_WT_M <- read_csv(
  "data/processed/DESeq2_KO_M_vs_WT_M_shrunk.csv",
  show_col_types = FALSE
)

HT_F_vs_WT_F <- read_csv(
  "data/processed/DESeq2_HT_F_vs_WT_F_shrunk.csv",
  show_col_types = FALSE
)

KO_M_vs_HT_F <- read_csv(
  "data/processed/DESeq2_KO_M_vs_HT_F_shrunk.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------------
# 4. Prepare volcano data
# ------------------------------------------------------------------

prepare_volcano_data <- function(deseq_df, comparison_name) {
  
  imprinted_genes_clean <- imprinted_genes %>%
    mutate(
      Gene_clean = str_to_upper(str_trim(Gene))
    )
  
  volcano_df <- deseq_df %>%
    mutate(
      Gene_Name = str_trim(Gene_Name),
      Gene_Base = str_remove(Gene_Name, "_129$|_JF1$"),
      Gene_clean = str_to_upper(Gene_Base)
    ) %>%
    left_join(
      imprinted_genes_clean,
      by = "Gene_clean"
    ) %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      padj > 0
    ) %>%
    mutate(
      Comparison = comparison_name,
      neg_log10_padj = -log10(padj),
      
      Is_Imprinted = !is.na(Expressed_Allele),
      Is_Significant = padj < PADJ_CUTOFF,
      Large_Fold_Change = abs(log2FoldChange) > LFC_CUTOFF,
      Significant_DEG = Is_Significant & Large_Fold_Change,
      Significant_Imprinted_DEG = Significant_DEG & Is_Imprinted,
      
      Volcano_Color = case_when(
        Significant_Imprinted_DEG & Expressed_Allele == "Maternal" ~ "Significant imprinted: maternal",
        Significant_Imprinted_DEG & Expressed_Allele == "Paternal" ~ "Significant imprinted: paternal",
        Significant_DEG ~ "Significant DEG",
        Is_Imprinted ~ "Imprinted but not significant",
        abs(log2FoldChange) <= LFC_CUTOFF ~ "Small fold change",
        TRUE ~ "Not significant"
      ),
      
      Label = case_when(
        Significant_Imprinted_DEG ~ Gene_Base,
        TRUE ~ NA_character_
      )
    )
  
  return(volcano_df)
}

# ------------------------------------------------------------------
# 5. Prepare volcano datasets
# ------------------------------------------------------------------

volcano_KO_M_vs_WT_M <- prepare_volcano_data(
  KO_M_vs_WT_M,
  "KO_M_vs_WT_M"
)

volcano_HT_F_vs_WT_F <- prepare_volcano_data(
  HT_F_vs_WT_F,
  "HT_F_vs_WT_F"
)

volcano_KO_M_vs_HT_F <- prepare_volcano_data(
  KO_M_vs_HT_F,
  "KO_M_vs_HT_F"
)

# ------------------------------------------------------------------
# 6. Make volcano plot function
# ------------------------------------------------------------------

make_volcano_plot <- function(volcano_df, title_text, output_file) {
  
  p <- ggplot(
    volcano_df,
    aes(x = log2FoldChange, y = neg_log10_padj)
  ) +
    geom_point(
      aes(color = Volcano_Color),
      alpha = 0.7,
      size = 1.6
    ) +
    
    geom_vline(
      xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.4
    ) +
    
    geom_hline(
      yintercept = -log10(PADJ_CUTOFF),
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.4
    ) +
    
    geom_text_repel(
      data = volcano_df %>% filter(!is.na(Label)),
      aes(label = Label, color = Volcano_Color),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50",
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Significant imprinted: maternal" = MATERNAL_COLOR,
        "Significant imprinted: paternal" = PATERNAL_COLOR,
        "Significant DEG" = "grey25",
        "Imprinted but not significant" = "darkgoldenrod2",
        "Small fold change" = "grey25",
        "Not significant" = "lightgray"
      )
    ) +
    
    labs(
      title = title_text,
      subtitle = paste0(
        "Significant DEGs: padj < ", PADJ_CUTOFF,
        " and |log2FC| > ", LFC_CUTOFF
      ),
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Gene category"
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    output_file,
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

# ------------------------------------------------------------------
# 7. Make volcano plots
# ------------------------------------------------------------------

volcano_plot_KO_M_vs_WT_M <- make_volcano_plot(
  volcano_KO_M_vs_WT_M,
  "Volcano Plot: KO Male vs WT Male",
  "figures/Volcano_KO_M_vs_WT_M.pdf"
)

volcano_plot_HT_F_vs_WT_F <- make_volcano_plot(
  volcano_HT_F_vs_WT_F,
  "Volcano Plot: HT Female vs WT Female",
  "figures/Volcano_HT_F_vs_WT_F.pdf"
)

volcano_plot_KO_M_vs_HT_F <- make_volcano_plot(
  volcano_KO_M_vs_HT_F,
  "Volcano Plot: KO Male vs HT Female",
  "figures/Volcano_KO_M_vs_HT_F.pdf"
)

# ------------------------------------------------------------------
# 8. Save significant imprinted DEG lists
# ------------------------------------------------------------------

sig_imprinted_KO_M_vs_WT_M <- volcano_KO_M_vs_WT_M %>%
  filter(Significant_Imprinted_DEG) %>%
  arrange(padj)

sig_imprinted_HT_F_vs_WT_F <- volcano_HT_F_vs_WT_F %>%
  filter(Significant_Imprinted_DEG) %>%
  arrange(padj)

sig_imprinted_KO_M_vs_HT_F <- volcano_KO_M_vs_HT_F %>%
  filter(Significant_Imprinted_DEG) %>%
  arrange(padj)

write_csv(
  sig_imprinted_KO_M_vs_WT_M,
  "data/processed/Sig_Imprinted_DEGs_KO_M_vs_WT_M.csv"
)

write_csv(
  sig_imprinted_HT_F_vs_WT_F,
  "data/processed/Sig_Imprinted_DEGs_HT_F_vs_WT_F.csv"
)

write_csv(
  sig_imprinted_KO_M_vs_HT_F,
  "data/processed/Sig_Imprinted_DEGs_KO_M_vs_HT_F.csv"
)

# ------------------------------------------------------------------
# 9. Read raw allele-specific counts and calculate AER
# ------------------------------------------------------------------

raw_counts <- read_csv(
  "data/raw/HIP_RNA-AlSp_rawcounts.csv",
  show_col_types = FALSE
)

sample_metadata <- tibble(
  Sample = c(
    "WT_M1", "WT_M2", "WT_M3",
    "WT_F1", "WT_F2", "WT_F3",
    "HT_F1", "HT_F2", "HT_F3",
    "KO_M1", "KO_M2", "KO_M3"
  ),
  Genotype = c(
    "WT_M", "WT_M", "WT_M",
    "WT_F", "WT_F", "WT_F",
    "HT_F", "HT_F", "HT_F",
    "KO_M", "KO_M", "KO_M"
  )
)

aer_long <- raw_counts %>%
  mutate(
    Gene_Base = str_remove(Gene_Name, "_129$|_JF1$"),
    Allele = case_when(
      str_detect(Gene_Name, "_129$") ~ "129",
      str_detect(Gene_Name, "_JF1$") ~ "JF1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Allele)) %>%
  pivot_longer(
    cols = all_of(names(sample_name_map)),
    names_to = "Original_Sample",
    values_to = "Raw_Count"
  ) %>%
  mutate(
    Raw_Count = as.numeric(Raw_Count),
    Sample = sample_name_map[Original_Sample],
    Allele_Type = case_when(
      Allele == MATERNAL_ALLELE ~ "Maternal",
      Allele == PATERNAL_ALLELE ~ "Paternal",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Allele_Type)) %>%
  group_by(Gene_Base, Sample, Original_Sample, Allele_Type) %>%
  summarise(
    Raw_Count = sum(Raw_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Allele_Type,
    values_from = Raw_Count
  ) %>%
  mutate(
    AER = case_when(
      Maternal + Paternal > 0 ~ (Maternal - Paternal) / (Maternal + Paternal),
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(sample_metadata, by = "Sample")

# ------------------------------------------------------------------
# 10. Make AER violin plot function
# ------------------------------------------------------------------

make_aer_violin_plot <- function(sig_gene_df, comparison_groups, title_text, output_file) {
  
  genes_to_plot <- sig_gene_df %>%
    mutate(
      Gene_Base = str_remove(Gene_Name, "_129$|_JF1$")
    ) %>%
    pull(Gene_Base) %>%
    unique()
  
  if (length(genes_to_plot) == 0) {
    message("No significant genes found for: ", title_text)
    return(NULL)
  }
  
  plot_data <- aer_long %>%
    filter(
      Gene_Base %in% genes_to_plot,
      Genotype %in% comparison_groups,
      !is.na(AER)
    ) %>%
    mutate(
      Genotype = factor(Genotype, levels = comparison_groups),
      Gene_Base = factor(Gene_Base, levels = genes_to_plot)
    )
  
  p <- ggplot(
    plot_data,
    aes(x = Genotype, y = AER, fill = Genotype)
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.45,
      color = "black",
      linewidth = 0.35
    ) +
    
    geom_jitter(
      aes(color = Genotype),
      width = 0.08,
      height = 0,
      size = 2.2,
      alpha = 0.95
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.4
    ) +
    
    geom_hline(
      yintercept = c(-0.7, 0.7),
      linetype = "dotted",
      color = "grey50",
      linewidth = 0.4
    ) +
    
    scale_y_continuous(
      limits = c(-1, 1),
      breaks = c(-1, -0.7, 0, 0.7, 1)
    ) +
    
    scale_fill_manual(
      values = setNames(c("#00274C", "#FFCB05"), comparison_groups)
    ) +
    
    scale_color_manual(
      values = setNames(c("black", "black"), comparison_groups)
    ) +
    
    facet_wrap(
      ~ Gene_Base,
      scales = "fixed",
      ncol = 4
    ) +
    
    labs(
      title = title_text,
      subtitle = "AER = (Maternal - Paternal) / (Maternal + Paternal)",
      x = "Genotype",
      y = "Allelic Expression Ratio (AER)"
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold", size = 8),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  print(p)
  
  ggsave(
    output_file,
    p,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  return(p)
}

# ------------------------------------------------------------------
# 11. Make violin plots
# ------------------------------------------------------------------

violin_KO_M_vs_WT_M <- make_aer_violin_plot(
  sig_imprinted_KO_M_vs_WT_M,
  comparison_groups = c("WT_M", "KO_M"),
  title_text = "Imprinted Significant DEGs: KO Male vs WT Male",
  output_file = "figures/Violin_Imprinted_DEGs_KO_M_vs_WT_M.jpg"
)

violin_HT_F_vs_WT_F <- make_aer_violin_plot(
  sig_imprinted_HT_F_vs_WT_F,
  comparison_groups = c("WT_F", "HT_F"),
  title_text = "Imprinted Significant DEGs: HT Female vs WT Female",
  output_file = "figures/Violin_Imprinted_DEGs_HT_F_vs_WT_F.jpg"
)

violin_KO_M_vs_HT_F <- make_aer_violin_plot(
  sig_imprinted_KO_M_vs_HT_F,
  comparison_groups = c("HT_F", "KO_M"),
  title_text = "Imprinted Significant DEGs: KO Male vs HT Female",
  output_file = "figures/Violin_Imprinted_DEGs_KO_M_vs_HT_F.jpg"
)

###############################################################################
# End of script
###############################################################################
