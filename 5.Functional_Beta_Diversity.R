
rm(list = ls())
setwd("the_path")
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
data_file <- "functional_data.xlsx"  # Must be in working directory (English name!)

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nPlease rename '功能数据库.xlsx' to 'functional_data.xlsx'.")
}

group_levels <- c("Control", "M7", "M14")

# =============================================================================
# Load and Prepare Data
# =============================================================================

# Load KEGG L3 pathway abundance
kegg_sheet <- "kegg_pathway_summary_L3"
kegg_raw <- readxl::read_xlsx(data_file, sheet = kegg_sheet)
rownames(kegg_raw) <- kegg_raw$Kegg3
kegg_t <- as.data.frame(t(kegg_raw[, -c(1, 2)]))  # Remove first two columns, transpose

# Load group metadata
group_raw <- readxl::read_xlsx(data_file, sheet = "Group")
rownames(group_raw) <- group_raw$ID

# Merge by sample ID (row names)
merged <- merge(group_raw, kegg_t, by = "row.names")
merged <- merged[merged$Group != "M28", ]  # Exclude M28 group
merged$Group <- factor(merged$Group, levels = group_levels)

# Final data: rows = samples, columns = pathways (first column is Group)
data_matrix <- merged[, -1]  # Numeric matrix
group_vector <- merged$Group

message("Sample counts per group:")
print(table(group_vector))

# =============================================================================
# 1. PCA using ropls
# =============================================================================
set.seed(123)
pca_model <- opls(x = data_matrix, predI = 2)

pca_scores <- as.data.frame(pca_model@scoreMN[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Group <- group_vector

var_pc1 <- round(100 * pca_model@modelDF[1, 1], 1)
var_pc2 <- round(100 * pca_model@modelDF[2, 1], 1)

# ⚠️ FIX: Y-axis was incorrectly labeled as "PC1" — now corrected to "PC2"
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95) +
  scale_color_d3() +
  labs(
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)"),  # ← FIXED!
    title = "PCA"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_pubr(legend = "bottom") +
  theme(plot.title = element_text(hjust = 0))

ggsave(p_pca, "p_pca.png", width = 6, height = 6, dpi = 300)
ggsave(p_pca, "p_pca.pdf", width = 6, height = 6)

# Save PCA coordinates
write.csv(pca_scores, "pca_coordinates.csv", row.names = FALSE)

# =============================================================================
# 2. PCoA (Bray-Curtis)
# =============================================================================
dist_bray <- vegdist(data_matrix, method = "bray")
pcoa_result <- cmdscale(dist_bray, eig = TRUE, k = 2)

# Calculate variance explained by positive eigenvalues only
eig_all <- pcoa_result$eig
eig_positive <- eig_all[eig_all > 0]
total_var <- sum(eig_positive)
pcoa_var <- round(100 * pcoa_result$eig[1:2] / total_var, 2)

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


