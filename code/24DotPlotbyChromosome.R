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
# 16. Manually assign gene order and genotype row order
# ------------------------------------------------------------------

manual_gene_order <- tibble(
  Gene = c(
    "Zdbf2",
    "Sfmbt2",
    "Wt1",
    "Gatm",
    "H13",
    "Mcts2",
    "Nnat",
    "Blcap",
    "Nespas",
    "Gnas",
    "Calcr",
    "Tfpi2",
    "Sgce",
    "Peg10",
    "Ppp1r9a",
    "Pon3",
    "Pon2",
    "Asb4",
    "Mest",
    "Copg2",
    "Klf14",
    "Nap1l5",
    "Zim2",
    "Zim1",
    "Peg3",
    "Usp29",
    "Zim3",
    "Zfp264",
    "Axl",
    "Ube3a",
    "Snrpn",
    "Ndn",
    "Magel2",
    "Mkrn3",
    "Zfp127as",
    "Peg12",
    "Ampd3",
    "Inpp5f",
    "H19",
    "Igf2as",
    "Igf2",
    "Ins2",
    "Th",
    "Ascl2",
    "Cd81",
    "Tssc4",
    "Kcnq1",
    "Kcnq1ot1",
    "Cdkn1c",
    "AF313042",
    "Slc22a18",
    "Phlda2",
    "Nap1l4",
    "Tnfrsf23",
    "Dhcr7",
    "Mir184",
    "Rasgrf1",
    "Plagl1",
    "Dcn",
    "Ddc",
    "Grb10",
    "Commd1",
    "Begain",
    "Dlk1",
    "Meg3",
    "Rtl1",
    "Rian",
    "Mirg",
    "Dio3",
    "Htr2a",
    "Kcnk9",
    "Peg13",
    "Trappc9",
    "Slc38a4",
    "Slc22a3",
    "Slc22a2",
    "Igf2r",
    "Impact"
    # add the rest manually in the order you want
  ),
  Gene_Order = c(
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78
  )
)

genotype_order <- c("WT Male", "KO Male", "WT Female", "HT Female")

row_gap <- 0.65

aer_plot_data <- aer_long %>%
  filter(!is.na(AER)) %>%
  inner_join(manual_gene_order, by = "Gene") %>%
  mutate(
    Group = factor(Group, levels = genotype_order),
    Group_Order = as.numeric(Group),
    
    # 4 rows per gene
    # Gene 1 = top, Gene 78 = bottom
    y_position = -((Gene_Order - 1) * (4 + row_gap) + Group_Order),
    
    y_label = paste0(Gene_Order, ". ", Gene, " — ", Group)
  )

gene_backgrounds <- aer_plot_data %>%
  group_by(Gene_Order, Gene) %>%
  summarize(
    ymin = min(y_position) - 0.5,
    ymax = max(y_position) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    background_fill = if_else(Gene_Order %% 2 == 1, "gray95", "white")
  )
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
  aes(x = AER, y = y_position)
) +
  
  geom_rect(
    data = gene_backgrounds,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax,
      fill = background_fill
    ),
    inherit.aes = FALSE,
    alpha = 1
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
  
  scale_fill_identity() +
  
  geom_point(
    aes(
      shape = Group,
      color = Observed_Expression,
      size = log2_total_count
    ),
    alpha = 0.85,
    position = position_jitter(width = 0.03, height = 0.06)
  ) +
  
  scale_shape_manual(
    values = c(
      "WT Male" = 16,
      "KO Male" = 17,
      "WT Female" = 15,
      "HT Female" = 18
    )
  ) +
  
  scale_y_continuous(
    breaks = aer_plot_data %>%
      distinct(y_position, y_label) %>%
      arrange(y_position) %>%
      pull(y_position),
    labels = aer_plot_data %>%
      distinct(y_position, y_label) %>%
      arrange(y_position) %>%
      pull(y_label)
  ) +
  
  scale_x_continuous(
    breaks = c(-1, -0.7, 0, 0.7, 1),
    labels = c("Paternal\n-1", "-0.7", "Balanced\n0", "0.7", "Maternal\n+1")
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  
  scale_color_manual(
    values = c(
      "Maternal" = "#FF6090",
      "Paternal" = "#1E4E79",
      "Not Imprinted" = "#B0B0B0",
      "No Expression" = "#E6E6E6"
    ),
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name = "log2(total count + 1)",
    range = c(1.8, 5)
  ) +
  
  labs(
    title = "Allele-Specific Expression Bias for Published Imprinted Genes",
    subtitle = "Each gene is manually ordered; genotype rows are WT Male, KO Male, WT Female, HT Female",
    x = "Allelic Expression Ratio",
    y = "Gene and genotype",
    color = "Observed Expression",
    shape = "Genotype"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

aer_dot_plot


# ------------------------------------------------------------------
# 20. Save AER dot plot
# ------------------------------------------------------------------
ggsave(
  "figures/Imprinted_Genes_AER_DotPlot_Chromosome_Ordered.jpg",
  aer_dot_plot,
  width = 18,
  height = 24,
  dpi = 300
)

##############################################################################
# End of script
##############################################################################