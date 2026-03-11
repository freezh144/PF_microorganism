rm(list = ls())

# Load required packages
required_pkgs <- c(
  "readxl", "tidyverse", "microeco", "ggplot2", "ggpubr", 
  "ggsci", "patchwork", "vegan"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================
data_file <- "microbiome_diff_analysis.xlsx"  # Must be in working directory

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nPlease rename your file to 'microbiome_diff_analysis.xlsx' or update 'data_file'.")
}

group_levels <- c("Control", "M7", "M14")

# =============================================================================
# Load Data Sheets
# =============================================================================

# OTU table (samples as rows, taxa as columns)
otu <- readxl::read_xlsx(data_file, sheet = "otu")
rownames(otu) <- otu$ID
otu <- otu[, -1]  # Remove ID column

# Group metadata
group <- readxl::read_xlsx(data_file, sheet = "group")
rownames(group) <- group$sample

# Taxonomy table (same row order as OTU; taxa as rows, ranks as columns)
tax <- readxl::read_xlsx(data_file, sheet = "tax")
rownames(tax) <- tax$ID
tax <- tax[, -1]

message("Sample counts per group:")
print(table(group$Group))

# =============================================================================
# Create microtable Object
# =============================================================================
dataset <- microtable$new(
  sample_table = group,
  otu_table = otu,
  tax_table = tax
)

message("microtable object created successfully.")

# =============================================================================
# Run LEfSe Analysis
# =============================================================================
lefse <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "Group",
  alpha = 0.3,                # Kruskal-Wallis test p-value threshold
  p_adjust_method = "fdr",    # Multiple testing correction
  lefse_subgroup = NULL
)

# Save results
write.csv(lefse$res_diff, "lefse_results.csv", row.names = FALSE)
message("LEfSe results saved to 'lefse_results.csv'.")

# =============================================================================
# Visualization 1: Top 30 Features by LDA Score (Bar Plot)
# =============================================================================
p_lefse_bar <- lefse$plot_diff_bar(
  use_number = 1:30,
  width = 0.8,
  group_order = group_levels
) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_d3() +
  scale_fill_d3() +
  theme_pubr(legend = "right")

ggsave(p_lefse_bar, "p_lefse_bar.png", width = 9, height = 6, dpi = 300)
ggsave(p_lefse_bar, "p_lefse_bar.pdf", width = 9, height = 6)

# =============================================================================
# Visualization 2: Cladogram (Top 200 Taxa, Top 50 Features)
# =============================================================================
p_lefse_cladogram <- lefse$plot_diff_cladogram(
  use_taxa_num = 200,
  use_feature_num = 50,
  clade_label_level = 5,
  group_order = group_levels
)

ggsave(p_lefse_cladogram, "p_lefse_cladogram.png", width = 8, height = 6, dpi = 300)
ggsave(p_lefse_cladogram, "p_lefse_cladogram.pdf", width = 8, height = 6)

