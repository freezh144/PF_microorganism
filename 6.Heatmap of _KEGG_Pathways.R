rm(list = ls())

# Load required packages
required_pkgs <- c("pheatmap", "ggplotify", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Load Data
# =============================================================================

# Expression matrix: rows = pathways, columns = samples
exp_file <- "dif_kegg_exp.csv"
if (!file.exists(exp_file)) stop("File not found: ", exp_file)
exp <- read.csv(exp_file, header = TRUE, row.names = 1)

# Group metadata: rows = samples, column = Group
group_file <- "Group.csv"
if (!file.exists(group_file)) stop("File not found: ", group_file)
Group <- read.csv(group_file, header = TRUE, row.names = 1)

# Ensure sample order matches between exp and Group
common_samples <- intersect(colnames(exp), rownames(Group))
if (length(common_samples) == 0) {
  stop("No matching sample IDs between expression matrix and group file.")
}
exp <- exp[, common_samples, drop = FALSE]
Group <- Group[common_samples, , drop = FALSE]

message("Sample counts per group:")
print(table(Group$Group))

# =============================================================================
# Define Colors and Annotations
# =============================================================================

# Custom diverging color palette: blue → white → orange
mycol <- colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200)

# Annotation colors for groups
# ⚠️ Ensure group labels in data EXACTLY match these names
ann_colors <- list(
  Group = c(
    Control = "#1F77B4FF",
    M07     = "#FF7F0EFF",
    M14     = "#2CA02CFF"
  )
)

# Convert Group to factor with explicit levels for consistent ordering
Group$Group <- factor(Group$Group, levels = c("Control", "M07", "M14"))

# =============================================================================
# Generate Heatmap
# =============================================================================

# Create pheatmap object
pheatmap_obj <- pheatmap(
  mat = as.matrix(exp),
  scale = "row",                     # Z-score per pathway
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = 3,                   # Optional: highlight 3 row clusters
  treeheight_row = 0,                # Hide dendrogram
  border_color = "white",
  color = mycol,
  annotation_col = Group,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 9,
  fontsize_number = 14,
  cellwidth = 15,
  cellheight = 15
)

# Convert to ggplot for saving via ggsave
library(ggplotify)
p_heatmap <- as.ggplot(pheatmap_obj)

# Save high-resolution figures
ggsave(p_heatmap, "pheatmap_KEGG.png", width = 12, height = 8, dpi = 300)
ggsave(p_heatmap, "pheatmap_KEGG.pdf", width = 12, height = 8)


