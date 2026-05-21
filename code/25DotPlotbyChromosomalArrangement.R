#############################################################################
# Author: Kylee Duczyminski
# About: This script creates a dot plot that arrangements the genes in order
# of their location on the mouse chromosomes. This allows us to see if there
# are clear differences in imprinting patterns at imprinted gene clusters.
#
# Input files:
#   data/raw/HIP_RNA-AlSp_rawcounts.csv
#   data/raw/imprinted_gene_list.csv
#
# Output files:
#   data/processed/HIP_RNA-AlSp_rawcounts_imprinted_AER_AllImprintedGenes.csv
#   figures/ALL_Imprinted_Genes_AER_DotPlot.jpg
#############################################################################

# Load libraries
library(tidyverse)

# Set working directory
setwd("/Users/kyleeduczyminski/Documents/Iwase-Lab/Kylee-Project")

# Read input files
raw_counts <- read_csv(
  "data/raw/HIP_RNA-AlSp_rawcounts.csv",
  show_col_types = FALSE
)

imprinted_genes <- read_csv(
  "data/raw/imprinted_gene_list.csv",
  show_col_types = FALSE
)

# Rename raw sequencing columns to biological sample names
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

# Define sample groups
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

# Clean published imprinted gene list

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

# Rename first column in raw counts
# First column contains values like:
#   Gene_129
#   Gene_JF1
names(raw_counts)[1] <- "Gene_Allele"

# Parse gene name and allele source
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

# Check whether parsing worked
cat("\nAllele source counts after parsing:\n")

raw_counts_parsed %>%
  count(Allele_Source) %>%
  print()

cat("\nFirst 10 parsed gene rows:\n")

raw_counts_parsed %>%
  select(Gene_Allele, Gene, Allele_Source) %>%
  head(10) %>%
  print()

# Check gene overlap between raw counts and published imprinted list
raw_gene_names <- unique(raw_counts_parsed$Gene)
published_gene_names <- unique(imprinted_genes$Gene)

overlapping_genes <- intersect(raw_gene_names, published_gene_names)
missing_from_raw <- setdiff(published_gene_names, raw_gene_names)

cat("\nNumber of unique genes in raw counts:", length(raw_gene_names), "\n")
cat("Number of genes in published imprinted list:", length(published_gene_names), "\n")
cat("Number of overlapping genes:", length(overlapping_genes), "\n")

cat("\nPublished imprinted genes missing from raw counts:\n")
print(missing_from_raw)

# Keep only imprinted published genes
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

# Convert to long format
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

# Convert maternal and paternal rows into columns
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

# Calculate total counts and AER
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

# Save AER long-format file
write_csv(
  aer_long,
  "data/processed/HIP_RNA-AlSp_rawcounts_imprinted_AER_long.csv"
)

# Manually assign gene order and genotype row order
manual_gene_order <- tibble::tibble(
  Gene = c(
    
    # Chr 1
    "Zdbf2",
    "Adam23",
    
    # Chr 2
    "Mcts2",
    "Mir296",
    "Sfmbt2",
    "Wt1",
    "Gatm",
    "H13",
    "Blcap",
    "Nnat",
    "Zfp64",
    "Gnas",
    
    # Chr 3
    "Jade1",
    "Gnai3",
    
    # Chr 5
    "Magi2",
    "Fkbp6",
    
    # Chr 6
    "Peg10",
    "Pon3",
    "Ppp1r9a",
    "Tfpi2",
    "Sgce",
    "Asb4",
    "Calcr",
    "Mest",
    "Aqp1",
    "Pon2",
    "Copg2",
    "Klf14",
    "Nap1l5",
    
    # Chr 7
    "Ano1",
    "Nlrp2",
    "Peg3os",
    "Usp29",
    "Peg3",
    "Zim1",
    "Zim3",
    "Axl",
    "Peg12",
    "Magel2",
    "Ndn",
    "Ube3a",
    "Mkrn3",
    "Snrpn",
    "Ampd3",
    "H19",
    "Ins2",
    "Kcnq1",
    "Ascl2",
    "Tssc4",
    "Cd81",
    "Kcnq1ot1",
    "Phlda2",
    "Slc22a18",
    "Cdkn1c",
    "Nap1l4",
    "Tnfrsf23",
    "Tnfrsf22",
    "Tnfrsf26",
    "Inpp5f",
    "Nctc1",
    "Igf2",
    "Th",
    "Dhcr7",
    "Zim2",
    "Snurf",
    
    # Chr 8
    "Gab1",
    "Sall1",
    "Cdh15",
    
    # Chr 9
    "Snx14",
    "Mir184",
    "Rasgrf1",
    
    # Chr 10
    "Hymai",
    "Plagl1",
    "Dcn",
    
    # Chr 11
    "Platr20",
    "Ddc",
    "Grb10",
    "Commd1",
    "Ccdc40",
    
    # Chr 12
    "Smoc1",
    "Dlk1",
    "Meg3",
    "Rian",
    "Dio3",
    "Mirg",
    "Rtl1",
    "Begain",
    
    # Chr 14
    "Htr2a",
    
    # Chr 15
    "Trappc9",
    "Chrac1",
    "Galnt6",
    "Peg13",
    "Slc38a4",
    "Kcnk9",
    
    # Chr 17
    "Qpct",
    "Pde10a",
    "Dact2",
    "Prkn",
    "Slc22a3",
    "Slc22a2",
    "Igf2r",
    "Thbs2",
    "Smoc2",
    
    # Chr 18
    "Impact",
    
    # Chr X
    "Fthl17b",
    "Fthl17c",
    "Fthl17d",
    "Gm35612",
    "Fthl17e",
    "Fthl17f",
    "Fthl17a",
    "Xlr3b",
    "Tsix",
    "Jpx",
    "Ftx",
    "Maged2",
    "Xist"
  ),
  
  Chromosome = c(
    rep("Chr 1", 2),
    rep("Chr 2", 10),
    rep("Chr 3", 2),
    rep("Chr 5", 2),
    rep("Chr 6", 13),
    rep("Chr 7", 36),
    rep("Chr 8", 3),
    rep("Chr 9", 3),
    rep("Chr 10", 3),
    rep("Chr 11", 5),
    rep("Chr 12", 8),
    rep("Chr 14", 1),
    rep("Chr 15", 6),
    rep("Chr 17", 9),
    rep("Chr 18", 1),
    rep("Chr X", 13)
  )
) %>%
  dplyr::mutate(Gene_Order = dplyr::row_number())

genotype_order <- c("WT Male", "KO Male", "WT Female", "HT Female")

row_gap <- 5
within_gene_gap <- 5

aer_plot_data <- aer_long %>%
  inner_join(manual_gene_order, by = "Gene") %>%
  mutate(
    Group = factor(Group, levels = genotype_order),
    Group_Order = as.numeric(Group),
    
    AER_for_plot = if_else(is.na(AER), 0, AER),
    
    y_position = -(
      (Gene_Order - 1) * ((length(genotype_order) - 1) * within_gene_gap + row_gap) +
        (Group_Order - 1) * within_gene_gap),
    
    y_label = paste0(Gene_Order, ". ", Gene, " — ", Group),
    
    Point_Color_Group = case_when(
      is.na(AER) ~ "No Expression",
      AER > -0.7 & AER < 0.7 ~ "Middle / Not Imprinted",
      TRUE ~ as.character(Expressed_Allele)
    )
  )

size_limits <- range(aer_plot_data$log2_total_count, na.rm = TRUE)

gene_backgrounds <- aer_plot_data %>%
  group_by(Gene_Order, Gene) %>%
  summarize(
    ymin = min(y_position) - 0.45,
    ymax = max(y_position) + 0.45,
    .groups = "drop"
  ) %>%
  mutate(
    background_fill = if_else(Gene_Order %% 2 == 1, "gray95", "white")
  )

# Calculate group means for optional summary dots
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

# Create AER dot plot
# Interpretation:
#   Left side  = paternal / JF1 expression
#   Right side = maternal / 129 expression
#
# Dot color = Expected Expression for the Gene:
#   Dark Pink  = Maternal
#   Dark Blue = Paternal
#   Gray = Not Imprinted
#   Light gray = No Expression (not shown on AER axis because AER is NA)
#
# Dot size:
#   Larger dots = higher total expression

aer_dot_plot <- ggplot(
  aer_plot_data,
  aes(x = AER_for_plot, y = y_position)
) +
  
  # Alternating gray/white gene backgrounds
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
  
  # Observed paternal imprinting zone: AER -1 to -0.7
  annotate(
    "rect",
    xmin = -Inf,
    xmax = -0.7,
    ymin = -Inf,
    ymax = Inf,
    fill = "#B9F9FF",
    alpha = 0.25
  ) +
  
  annotate(
    "rect",
    xmin = 0.7,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "#FFD6EC",
    alpha = 0.25
  ) +
  
  # Center line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  
  # Imprinting threshold lines
  geom_vline(
    xintercept = c(-0.7, 0.7),
    linetype = "dotted",
    color = "gray30",
    linewidth = 0.5
  ) +
  
  scale_fill_identity() +
  
  # Dots are colored by EXPECTED expression from published list
  geom_point(
    aes(
      shape = Group,
      color = Point_Color_Group,
      size = log2_total_count
    ),
    alpha = 0.9,
    position = position_jitter(width = 0.03, height = 0.06)
  ) +
  
  scale_shape_manual(
    values = c(
      "WT Male" = 15,
      "KO Male" = 17,
      "WT Female" = 16,
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
    labels = c(
      "Paternal\n-1",
      "-0.7",
      "Balanced\n0",
      "0.7",
      "Maternal\n+1"
    )
  ) +
  
  coord_cartesian(xlim = c(-1, 1)) +
  
  # Expected expression colors
  scale_color_manual(
    name = "Dot Color",
    values = c(
      "Maternal" = "#C2185B",
      "Paternal" = "#0D47A1",
      "Middle / Not Imprinted" = "#B0B0B0",
      "No Expression" = "#E6E6E6"
    ),
    drop = FALSE
  ) +
  
  scale_size_continuous(
    name = "log2(total count + 1)",
    range = c(1.8, 5),
    limits = size_limits
  ) +
  
  labs(
    title = "Allele-Specific Expression Bias for Published Imprinted Genes",
    subtitle = "Dot color shows expected published expression; shaded regions show observed imprinting threshold zones",
    x = "Allelic Expression Ratio",
    y = "Gene and genotype",
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

# Save AER dot plot
ggsave(
  "figures/ALL_Imprinted_Genes_AER_DotPlot_Chromosome_Ordered.pdf",
  aer_dot_plot,
  width = 6,
  height = 60,
  dpi = 300,
  limitsize = FALSE
)

# ===========================================================================
# DOT PLOT PER CHROMOSOME
# ===========================================================================

make_chr_plot <- function(chr) {
  
  chr_data <- aer_plot_data %>%
    filter(Chromosome == chr)
  
  chr_backgrounds <- gene_backgrounds %>%
    left_join(
      manual_gene_order %>% select(Gene, Gene_Order, Chromosome),
      by = c("Gene", "Gene_Order")
    ) %>%
    filter(Chromosome == chr)
  
  ggplot(chr_data, aes(x = AER_for_plot, y = y_position)) +
    geom_rect(
      data = chr_backgrounds,
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
    annotate(
      "rect",
      xmin = -Inf,
      xmax = -0.7,
      ymin = -Inf,
      ymax = Inf,
      fill = "#B9F9FF",
      alpha = 0.25
    ) +
    annotate(
      "rect",
      xmin = 0.7,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = "#FFD6EC",
      alpha = 0.25
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.7, 0.7), linetype = "dotted", color = "gray30", linewidth = 0.5) +
    scale_fill_identity() +
    geom_point(
      aes(
        shape = Group,
        color = Point_Color_Group,
        size = log2_total_count
      ),
      alpha = 0.9,
      position = position_jitter(width = 0.03, height = 0.06)
    ) +
    scale_shape_manual(
      values = c(
        "WT Male" = 15,
        "KO Male" = 17,
        "WT Female" = 16,
        "HT Female" = 18
      )
    ) +
    scale_y_continuous(
      breaks = chr_data %>%
        distinct(y_position, y_label) %>%
        arrange(y_position) %>%
        pull(y_position),
      labels = chr_data %>%
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
      name = "Dot Color",
      values = c(
        "Maternal" = "#C2185B", ##C2185B
        "Paternal" = "#0D47A1", ##0D47A1
        "Middle / Not Imprinted" = "#B0B0B0",
        "No Expression" = "#E6E6E6"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(
      name = "log2(total count + 1)",
      range = c(1.8, 5),
      limits = size_limits
    ) +
    labs(
      title = paste0("Allele-Specific Expression Bias — ", chr),
      subtitle = "Threshold-crossing dots are colored by expected expression; middle/no-expression dots are gray",
      x = "Allelic Expression Ratio",
      y = "Gene and genotype",
      shape = "Genotype"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
}

height_per_gene <- 2
base_height <- 2

for (chr in unique(manual_gene_order$Chromosome)) {
  
  chr_plot <- make_chr_plot(chr)
  
  safe_chr <- str_replace_all(chr, " ", "_")
  
  output_file <- file.path(
    "figures",
    paste0("AER_DotPlot_", safe_chr, ".pdf")
  )
  
  n_genes_chr <- n_distinct(
    manual_gene_order$Gene[manual_gene_order$Chromosome == chr]
  )
  
  message("Saving: ", output_file)
  
  ggsave(
    filename = output_file,
    plot = chr_plot,
    width = 8,
    height = base_height + n_genes_chr * height_per_gene,
    device = "pdf",
    limitsize = FALSE
  )
}
