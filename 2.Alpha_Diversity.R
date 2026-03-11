# Clear workspace
rm(list = ls())

# Load required packages
required_pkgs <- c("ggplot2", "ggpubr", "readxl", "ggsci", "patchwork")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================
data_file <- "pulmonary_fibrosis_metagenomics.xlsx"  # Must be in working directory

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file, 
       "\nPlease ensure the file exists and is named in English.")
}

# =============================================================================
# Load and Prepare Data
# =============================================================================
# Read alpha diversity metrics (Sheet 4)
alpha_data <- readxl::read_xlsx(data_file, sheet = 4)
alpha_data <- as.data.frame(alpha_data)

# Read group metadata (Sheet named "Group")
group_data <- readxl::read_xlsx(data_file, sheet = "Group")
group_data <- as.data.frame(group_data)

# Merge by sample ID
# Assuming: group_data$ID matches alpha_data$sample
merged_data <- merge(
  x = group_data,
  y = alpha_data,
  by.x = "ID",
  by.y = "sample",
  all = FALSE
)

# Remove M28 group and set factor levels
merged_data <- subset(merged_data, Group != "M28")
merged_data$Group <- factor(merged_data$Group, levels = c("Control", "M7", "M14"))

message("Sample counts per group:")
print(table(merged_data$Group))

# =============================================================================
# Helper Function: Create Violin + Box Plot with Stats
# =============================================================================
create_alpha_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = Group, y = .data[[y_var]])) +
    geom_violin(aes(fill = Group, color = Group), alpha = 0.8) +
    geom_boxplot(
      aes(x = Group, y = .data[[y_var]]),
      fill = "white",
      color = "gray",
      width = 0.08,
      outlier.shape = NA
    ) +
    scale_fill_d3() +
    scale_color_d3() +
    labs(x = NULL, y = y_label, title = NULL) +
    theme_pubr(legend = "none") +
    ggpubr::stat_compare_means(
      method = "anova",  # global test
      label = "p.format",
      tip.length = 0.01
    )
}

# =============================================================================
# Generate Individual Plots
# =============================================================================
p_simpson <- create_alpha_plot(merged_data, "simpson", "Simpson")
p_chao1   <- create_alpha_plot(merged_data, "chao1",   "Chao1")
p_ace     <- create_alpha_plot(merged_data, "ACE",     "ACE")
p_shannon <- create_alpha_plot(merged_data, "shannon", "Shannon")

# =============================================================================
# Combine Plots in Different Layouts
# =============================================================================

# Layout 1: 2x2 grid
p_alpha_grid <- p_simpson + p_chao1 + p_ace + p_shannon +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")
ggsave(p_alpha_grid, "p_alpha_grid.png", width = 12, height = 12, dpi = 300)
ggsave(p_alpha_grid, "p_alpha_grid.pdf", width = 12, height = 12)

# Layout 2: Horizontal row (all four)
p_alpha_row <- p_simpson + p_chao1 + p_ace + p_shannon +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A")
ggsave(p_alpha_row, "p_alpha_row.png", width = 15, height = 5, dpi = 300)
ggsave(p_alpha_row, "p_alpha_row.pdf", width = 15, height = 5)

# Layout 3: Only Simpson and Shannon (commonly reported)
p_alpha_two <- p_simpson + p_shannon +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A")
ggsave(p_alpha_two, "p_alpha_two.png", width = 7, height = 5, dpi = 300)
ggsave(p_alpha_two, "p_alpha_two.pdf", width = 7, height = 5)


