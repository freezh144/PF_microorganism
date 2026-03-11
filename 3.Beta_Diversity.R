rm(list = ls())

# Load required packages
required_pkgs <- c("ggplot2", "ggpubr", "readxl", "ropls", "vegan", "patchwork", "ggsci")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================
data_file <- "taxa_data.xlsx"  # Must be in working directory (use English name!)

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nPlease rename your file to 'taxa_data.xlsx' or update 'data_file'.")
}

# Group levels for consistent ordering
group_levels <- c("Control", "M7", "M14")

# =============================================================================
# Load and Prepare Data
# =============================================================================

# Load species abundance (L7/Species level)
species_sheet <- "species_annotation_L7_Species"
abundance_raw <- readxl::read_xlsx(data_file, sheet = species_sheet)
rownames(abundance_raw) <- abundance_raw$Species
abundance_t <- as.data.frame(t(abundance_raw[, -c(1, 2)]))  # Remove first two columns, transpose

# Load group metadata
group_raw <- readxl::read_xlsx(data_file, sheet = "Group")
rownames(group_raw) <- group_raw$ID

# Merge by row names (sample IDs)
merged <- merge(group_raw, abundance_t, by = "row.names")
merged <- merged[merged$Group != "M28", ]  # Exclude M28
merged$Group <- factor(merged$Group, levels = group_levels)

# Final data matrix: rows = samples, columns = taxa (first column is Group)
data_matrix <- merged[, -1]  # Abundance matrix (numeric)
group_vector <- merged$Group  # Factor vector

message("Sample counts per group:")
print(table(group_vector))

# =============================================================================
# 1. PCA using ropls
# =============================================================================
set.seed(123)
pca_model <- opls(x = data_matrix, predI = 2)

# Extract scores and variance explained
pca_scores <- as.data.frame(pca_model@scoreMN[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Group <- group_vector

var_pc1 <- round(100 * pca_model@modelDF[1, 1], 1)
var_pc2 <- round(100 * pca_model@modelDF[2, 1], 1)

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95) +
  scale_color_d3() +
  labs(
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)"),
    title = "PCA"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_pubr(legend = "bottom") +
  theme(plot.title = element_text(hjust = 0))

ggsave(p_pca, "p_pca.png", width = 6, height = 6, dpi = 300)
ggsave(p_pca, "p_pca.pdf", width = 6, height = 6)

# =============================================================================
# 2. PCoA (Bray-Curtis)
# =============================================================================
dist_bray <- vegdist(data_matrix, method = "bray")
pcoa_result <- cmdscale(dist_bray, eig = TRUE, k = 2)

# Eigenvalues and variance explained
eig_vals <- pcoa_result$eig[1:2]
total_var <- sum(pcoa_result$eig[pcoa_result$eig > 0])
pcoa_var <- round(100 * eig_vals / total_var, 2)

pcoa_df <- data.frame(
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2],
  Group = group_vector
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(fill = Group), geom = "polygon", level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_color_d3() +
  scale_fill_d3() +
  labs(
    x = paste0("PCoA1 (", pcoa_var[1], "%)"),
    y = paste0("PCoA2 (", pcoa_var[2], "%)"),
    title = "PCoA"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_pubr(legend = "right") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_pcoa, "p_pcoa.png", width = 7, height = 6, dpi = 300)
ggsave(p_pcoa, "p_pcoa.pdf", width = 7, height = 6)

# =============================================================================
# 3. NMDS (Bray-Curtis)
# =============================================================================
nmds_model <- metaMDS(dist_bray, k = 2, trymax = 100)
stress_val <- round(nmds_model$stress, 3)

nmds_df <- data.frame(
  NMDS1 = nmds_model$points[, 1],
  NMDS2 = nmds_model$points[, 2],
  Group = group_vector
)

p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group, shape = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(fill = Group), geom = "polygon", level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_color_d3() +
  scale_fill_d3() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = paste("NMDS\nStress =", stress_val)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_pubr(legend = "right") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_nmds, "p_nmds.png", width = 6, height = 6, dpi = 300)
ggsave(p_nmds, "p_nmds.pdf", width = 6, height = 6)

# =============================================================================
# Combine Plots
# =============================================================================

# Full panel: PCA + PCoA + NMDS
p_beta_full <- p_pca + p_pcoa + p_nmds +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(p_beta_full, "p_beta_full.png", width = 12, height = 5, dpi = 300)
ggsave(p_beta_full, "p_beta_full.pdf", width = 12, height = 5)

# Compact panel: PCoA + NMDS only
p_beta_compact <- p_pcoa + p_nmds +
  plot_layout(guides = "collect")
ggsave(p_beta_compact, "p_beta_compact.png", width = 10, height = 5, dpi = 300)
ggsave(p_beta_compact, "p_beta_compact.pdf", width = 10, height = 5)


