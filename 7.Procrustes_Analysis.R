
rm(list = ls())

# Set working directory to an ENGLISH PATH (e.g., F:/project/procrustes/)
setwd(getwd())

# Load required packages
required_pkgs <- c("vegan", "readxl", "ggplot2", "ggpubr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================
species_file <- "species_data.xlsx"      # Renamed from "物种数据表.xlsx"
function_file <- "functional_data.xlsx"  # Renamed from "功能数据库.xlsx"

if (!file.exists(species_file)) stop("Species file not found: ", species_file)
if (!file.exists(function_file)) stop("Functional file not found: ", function_file)

# =============================================================================
# Load and Prepare Data
# =============================================================================

# Load species data (rows = taxa, columns = samples)
species_raw <- readxl::read_xlsx(species_file, sheet = "species_annotation_L7_Species")
rownames(species_raw) <- species_raw$Species
species_t <- t(species_raw[, -c(1, 2)])  # Remove first two columns, transpose

# Load KEGG L3 pathway data (rows = pathways, columns = samples)
kegg_raw <- readxl::read_xlsx(function_file, sheet = "kegg_pathway_summary_L3")
rownames(kegg_raw) <- kegg_raw$Kegg3
kegg_t <- t(kegg_raw[, -c(1, 2)])

# Ensure sample IDs match between datasets
common_samples <- intersect(colnames(species_t), colnames(kegg_t))
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between species and functional data.")
}
species_t <- species_t[, common_samples]
kegg_t <- kegg_t[, common_samples]

message("Number of shared samples: ", ncol(species_t))

# =============================================================================
# Distance Calculation & NMDS
# =============================================================================

# Species: Bray-Curtis distance
spe_dist <- vegdist(species_t, method = "bray")

# Functions: Euclidean distance on scaled data
func_dist <- vegdist(scale(kegg_t), method = "euclidean")

# NMDS ordination (2D)
set.seed(123)
mds_species <- monoMDS(spe_dist, k = 2)
mds_function <- monoMDS(func_dist, k = 2)

# =============================================================================
# Procrustes Analysis & Significance Test
# =============================================================================

# Symmetric Procrustes rotation
procrustes_fit <- procrustes(mds_species, mds_function, symmetric = TRUE)

# PROTEST permutation test (999 permutations)
set.seed(1)
protest_result <- protest(mds_species, mds_function, permutations = 999)

m2_stat <- round(protest_result$ss, 3)
p_value <- protest_result$signif

message("Procrustes M² = ", m2_stat, ", p-value = ", p_value)

# =============================================================================
# Prepare Data for Plotting
# =============================================================================

# Extract coordinates
pro_data <- data.frame(
  X1 = procrustes_fit$Yrot[, 1],  # Rotated functional coords
  X2 = procrustes_fit$Yrot[, 2],
  MDS1 = procrustes_fit$X[, 1],   # Original species coords
  MDS2 = procrustes_fit$X[, 2]
)

rotation_matrix <- procrustes_fit$rotation

# =============================================================================
# Generate Procrustes Plot
# =============================================================================

p_procrustes <- ggplot(pro_data) +
  # Arrows: from species → midpoint → rotated function
  geom_segment(
    aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
    arrow = arrow(length = unit(0, 'cm')),
    color = "#9BBB59", size = 1
  ) +
  geom_segment(
    aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2),
    arrow = arrow(length = unit(0.2, 'cm')),
    color = "#957DB1", size = 1
  ) +
  # Points
  geom_point(aes(X1, X2), color = "#9BBB59", size = 3, shape = 16) +
  geom_point(aes(MDS1, MDS2), color = "#957DB1", size = 3, shape = 16) +
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # Rotation guide lines (optional)
  geom_abline(intercept = 0, slope = rotation_matrix[1, 2] / rotation_matrix[1, 1], 
              linetype = "dotted", color = "gray50") +
  geom_abline(intercept = 0, slope = rotation_matrix[2, 2] / rotation_matrix[2, 1], 
              linetype = "dotted", color = "gray50") +
  # Labels and title
  labs(
    x = "Dimension 1",
    y = "Dimension 2",
    title = "Procrustes Analysis: Species vs. KEGG Function"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf, hjust = 0, vjust = 1,
    label = paste0("M² = ", m2_stat, "\np = ", p_value),
    size = 4, fontface = "bold"
  ) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save figures
ggsave(p_procrustes, "procrustes_plot.png", width = 6, height = 6, dpi = 300)
ggsave(p_procrustes, "procrustes_plot.pdf", width = 6, height = 6)



