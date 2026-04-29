##############################################################################
# Author: Kylee Duczyminski
# About:
# This script creates TWO separate heatmaps:
#   1. Maternally expressed imprinted genes
#   2. Paternally expressed imprinted genes
#
# Each heatmap keeps the same sample order, color scheme, genotype annotation,
# and row clustering.
#
# Input:
# data/processed/HIP_RNA-AlSp_parentalASE_imprinting_analysis.csv
#
# Output:
# figures/Maternally_Expressed_Genes_Heatmap.pdf
# figures/Paternally_Expressed_Genes_Heatmap.pdf
##############################################################################

# ------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------
library(tidyverse)
library(pheatmap)

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
# 4. Define sample order
# ------------------------------------------------------------------
ordered_samples <- c(
  "WT_M1", "WT_M2", "WT_M3",
  "WT_F1", "WT_F2", "WT_F3",
  "HT_F1", "HT_F2", "HT_F3",
  "KO_M1", "KO_M2", "KO_M3"
)

# Build matching AER column names
ordered_aer_columns <- paste0(ordered_samples, "_AER")

# ------------------------------------------------------------------
# 5. Clean Expressed_Allele column
# ------------------------------------------------------------------
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
# 6. Keep only needed columns
# ------------------------------------------------------------------
imprinting_data <- imprinting_data %>%
  select(Gene_Name, Expressed_Allele, all_of(ordered_aer_columns))

# ------------------------------------------------------------------
# 7. Split into maternal and paternal datasets
# ------------------------------------------------------------------
maternal_data <- imprinting_data %>%
  filter(Expressed_Allele == "Maternal") %>%
  arrange(Gene_Name)

paternal_data <- imprinting_data %>%
  filter(Expressed_Allele == "Paternal") %>%
  arrange(Gene_Name)

# ------------------------------------------------------------------
# 8. Define heatmap function
# ------------------------------------------------------------------
make_imprinting_heatmap <- function(data, output_file, plot_title) {
  
  # Drop unused factor levels so the legend only shows the group being plotted
  data <- data %>%
    mutate(
      Expressed_Allele = droplevels(factor(Expressed_Allele))
    )
  
  # Row annotation: shows whether genes are maternally or paternally expressed
  annotation_row <- data %>%
    select(Gene_Name, Expressed_Allele) %>%
    column_to_rownames("Gene_Name")
  
  # Build numeric matrix of AER values
  mat <- data %>%
    select(Gene_Name, all_of(ordered_aer_columns)) %>%
    column_to_rownames("Gene_Name") %>%
    as.matrix()
  
  # Column annotation: sample genotype/sex group
  annotation_col <- data.frame(
    Genotype = c(
      rep("WT_M", 3),
      rep("WT_F", 3),
      rep("HT_F", 3),
      rep("KO_M", 3)
    )
  )
  rownames(annotation_col) <- colnames(mat)
  
  # Function to categorize AER values into color bins
  categorize_ratio <- function(x) {
    case_when(
      is.na(x) ~ "NA",
      x >= 0.7 ~ "red",
      x <= -0.7 ~ "blue",
      x > 0 & x < 0.7 ~ "light_red",
      x > -0.7 & x < 0 ~ "light_blue",
      x == 0 ~ "white"
    )
  }
  
  # Convert AER matrix into categorical matrix
  cat_mat <- matrix(
    categorize_ratio(as.vector(mat)),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  # Convert categories to numbers because pheatmap needs numeric values
  category_to_number <- c(
    "blue" = 1,
    "light_blue" = 2,
    "white" = 3,
    "light_red" = 4,
    "red" = 5,
    "NA" = 6
  )
  
  plot_mat <- matrix(
    category_to_number[cat_mat],
    nrow = nrow(cat_mat),
    ncol = ncol(cat_mat),
    dimnames = dimnames(cat_mat)
  )
  
  # Heatmap colors:
  # blue = strong paternal bias
  # light blue = slight paternal bias
  # white = no bias
  # light red = slight maternal bias
  # red = strong maternal bias
  # gray = NA
  colors <- c(
    "blue",
    "lightskyblue",
    "white",
    "lightcoral",
    "red",
    "gray"
  )
  
  breaks <- seq(0.5, 6.5, by = 1)
  
  # Annotation colors
  annotation_colors <- list(
    Genotype = c(
      WT_M = "#1f77b4",
      WT_F = "#ff7f0e",
      HT_F = "#2ca02c",
      KO_M = "#d62728"
    ),
    Expressed_Allele = c(
      Maternal = "#e41a1c",
      Paternal = "#377eb8"
    )
  )
  
  # Make heatmap
  pheatmap(
    plot_mat,
    color = colors,
    breaks = breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    border_color = "black",
    main = plot_title,
    filename = output_file,
    cellwidth = 18,
    cellheight = 8
  )
}

# ------------------------------------------------------------------
# 9. Create maternal heatmap
# ------------------------------------------------------------------
make_imprinting_heatmap(
  maternal_data,
  "figures/Maternally_Expressed_Genes_Heatmap.pdf",
  "Maternally Expressed Genes"
)

# ------------------------------------------------------------------
# 10. Create paternal heatmap
# ------------------------------------------------------------------
make_imprinting_heatmap(
  paternal_data,
  "figures/Paternally_Expressed_Genes_Heatmap.pdf",
  "Paternally Expressed Genes"
)



