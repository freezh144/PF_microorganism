# Plot Top 10 Taxa at Species, Genus, and Family Levels

# Clear workspace
rm(list = ls())

# Load required packages
required_packages <- c("ggplot2", "ggpubr", "reshape2", "ggalluvial", "patchwork", "readxl")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package '", pkg, "' is not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================
data_file <- "taxa_data_group.xlsx"
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

# Color palette (supports up to 25 taxa)
color_palette <- c(
  "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1",
  "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkgoldenrod1", "darkseagreen",
  "darkorchid", "#E64B35B2", "#4DBBD5B2", "#3C5488B2", "#3C5488B2", "#00A087B2",
  "#F39B7FB2", "#91D1C2B2", "#8491B4B2", "#DC0000B2", "#7E6148B2", "yellow", "darkolivegreen1"
)

# Plot dimensions
width_combined <- 12
height_combined <- 6
width_grouped <- 14
height_grouped <- 6

# =============================================================================
# Helper Function: Process and Plot Top 10 Taxa for a Given Level
# Arguments:
#   sheet_index: sheet number in Excel file (1=Species, 2=Genus, 3=Family)
#   taxon_level: character string ("Species", "Genus", or "Family")
# Returns:
#   List of ggplot objects: bar, flow, and combined
# =============================================================================
plot_top_taxa <- function(sheet_index, taxon_level) {
  # Read and preprocess data
  df_raw <- read_xlsx(data_file, sheet = sheet_index)
  df <- as.data.frame(df_raw)
  
  expected_cols <- c(taxon_level, "con_28D", "mod_7D", "mod_14D")
  if (!all(expected_cols %in% names(df))) {
    stop("Missing expected columns in sheet ", sheet_index)
  }
  
  df <- df[, expected_cols]
  names(df) <- c(taxon_level, "Control", "M7", "M14")
  df_top10 <- head(df, 10)
  
  # Reshape to long format
  df_long <- melt(
    df_top10,
    id.vars = taxon_level,
    variable.name = "Group",
    value.name = "Abundance"
  )
  
  # Stacked bar plot
  p_bar <- ggplot(df_long, aes(x = Group, y = 100 * Abundance, fill = .data[[taxon_level]])) +
    geom_col(position = "stack") +
    labs(x = NULL, y = "Relative abundance (%)", title = paste0("Top10_", taxon_level)) +
    scale_fill_manual(values = color_palette) +
    theme_pubr(legend = "right") +
    theme(
      legend.title = element_text(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black", size = 9),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(title = taxon_level, reverse = TRUE))
  
  # Alluvial plot
  p_flow <- ggplot(df_long, aes(
    x = Group,
    y = 100 * Abundance,
    fill = .data[[taxon_level]],
    stratum = .data[[taxon_level]],
    alluvium = .data[[taxon_level]]
  )) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "black", alpha = 0.4, width = 0.5) +
    geom_stratum(width = 0.5) +
    labs(x = NULL, y = "Relative abundance (%)", title = paste0("Top10_", taxon_level)) +
    scale_fill_manual(values = color_palette) +
    theme_pubr(legend = "right") +
    theme(
      legend.title = element_text(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black", size = 9),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(title = taxon_level, reverse = TRUE))
  
  # Combine
  p_combined <- p_bar + p_flow + plot_layout(guides = "collect")
  
  # Save outputs
  ggsave(p_combined, filename = paste0("p_", taxon_level, ".png"),
         width = width_combined, height = height_combined, dpi = 300)
  ggsave(p_combined, filename = paste0("p_", taxon_level, ".pdf"),
         width = width_combined, height = height_combined)
  
  return(list(bar = p_bar, flow = p_flow, combined = p_combined))
}

# =============================================================================
# Main Execution
# =============================================================================

# Generate plots for each taxonomic level
plots_species <- plot_top_taxa(sheet_index = 1, taxon_level = "Species")
plots_genus   <- plot_top_taxa(sheet_index = 2, taxon_level = "Genus")
plots_family  <- plot_top_taxa(sheet_index = 3, taxon_level = "Family")

# Group bar plots
grouped_bars <- plots_family$bar + plots_genus$bar + plots_species$bar +
  plot_annotation(tag_levels = "A")
ggsave(grouped_bars, "p_bar_plots.png", width = width_grouped, height = height_grouped, dpi = 300)
ggsave(grouped_bars, "p_bar_plots.pdf", width = width_grouped, height = height_grouped)

# Group alluvial plots
grouped_flows <- plots_family$flow + plots_genus$flow + plots_species$flow +
  plot_annotation(tag_levels = "A")
ggsave(grouped_flows, "p_alluvial_plots.png", width = width_grouped, height = height_grouped, dpi = 300)
ggsave(grouped_flows, "p_alluvial_plots.pdf", width = width_grouped, height = height_grouped)



