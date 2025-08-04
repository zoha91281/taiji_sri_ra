# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(grid)

# 1. Data preperation
# Load expression data and differential expression results
expr_df <- read.csv("rnaseq_updated_names.csv")
wilcox_df <- read.csv("wilcoxon_results_final.csv")

# Define sample groups
oa_samples <- c('OA1027', 'OA6', 'OA1', 'OA7', 'OA2', 'OA8', 'OA3', 'OA4', 'OA5', 'OA9', 'OA10')
ra_samples <- c('RA1', 'RA2', 'RA7', 'RA8', 'RA3', 'RA9', 'RA4', 'RA5', 'RA952', 'RA10', 'RA6')
all_samples <- c(oa_samples, ra_samples)

# Create annotation data frame
group_factor <- factor(ifelse(grepl("^OA", all_samples), "OA", "RA"))
annotation_col <- data.frame(Group = group_factor)
rownames(annotation_col) <- all_samples
palette <- c("OA" = "#74add1", "RA" = "#d73027")

# 2. Top 50 variable differntially expressed genes
# Filter significant genes (FDR < 0.05 and |log2FC| > log2(1.2))
significant <- wilcox_df[wilcox_df$FDR < 0.05 & abs(wilcox_df$log2FC_RA_vs_OA) > log2(1.2), ]
sig_genes <- significant$gene_name

# Extract expression data for significant genes
data_sig <- expr_df[expr_df$gene_name %in% sig_genes, c("gene_name", all_samples)]
rownames(data_sig) <- data_sig$gene_name
data_sig$gene_name <- NULL

# Select top 50 most variable genes
variances <- apply(data_sig, 1, var)
top_var_sig_genes <- names(sort(variances, decreasing = TRUE))[1:50]

# 3. Top 50 cytokine/TLR genes
# Identify cytokine and TLR genes (IL, TNF, TLR)
il_mask <- grepl("^IL", expr_df$gene_name, ignore.case = TRUE)
tnf_mask <- grepl("^TNF", expr_df$gene_name, ignore.case = TRUE)
tlr_mask <- grepl("^TLR", expr_df$gene_name, ignore.case = TRUE)
pathway_genes_all <- expr_df$gene_name[il_mask | tnf_mask | tlr_mask]

# Extract expression data
data_pathways <- expr_df[expr_df$gene_name %in% pathway_genes_all, c("gene_name", all_samples)]
rownames(data_pathways) <- data_pathways$gene_name
data_pathways$gene_name <- NULL

# Select top 50 most variable genes
pathway_variances <- apply(data_pathways, 1, var)
top_var_pathway_genes <- names(sort(pathway_variances, decreasing = TRUE))[1:50]

# 4. Heatmap Functions
# Z-score normalization function
zscore <- function(x) {
 (x - mean(x)) / sd(x)
}

# Heatmap plotting function
plot_heatmap <- function(data_matrix, annotation_col, title_text) {
 # Filter rows with variance > 0
 data_matrix <- data_matrix[apply(data_matrix, 1, var) > 0, ]
  # Z-score normalization and clipping
 data_scaled <- t(apply(data_matrix, 1, zscore))
 data_scaled[data_scaled > 2] <- 2
 data_scaled[data_scaled < -2] <- -2
  # Color palette
 heatmap_colors <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
  # Cluster columns separately for OA and RA
 gap_col_idx <- length(oa_samples)
 cluster_cols_oa <- hclust(as.dist(1 - cor(data_scaled[, oa_samples])))
 cluster_cols_ra <- hclust(as.dist(1 - cor(data_scaled[, ra_samples])))
 col_order <- c(oa_samples[cluster_cols_oa$order], ra_samples[cluster_cols_ra$order])
  # Create heatmap
 pheatmap(
   data_scaled[, col_order],
   color = heatmap_colors,
   cluster_rows = TRUE,
   cluster_cols = FALSE,
   show_rownames = TRUE,
   show_colnames = TRUE,
   annotation_col = annotation_col[col_order, , drop = FALSE],
   annotation_colors = list(Group = palette),
   border_color = NA,
   fontsize = 8,
   gaps_col = gap_col_idx,
   main = title_text
 )
}

# 5. Generate and Save Heatmaps
# Plot top 50 variable DEGs
plot_heatmap(
 data_sig[top_var_sig_genes, all_samples],
 annotation_col,
 "Top 50 Variable Differentially Expressed Genes (FDR < 0.05, |log2FC| > 1.2)"
)

# Plot top 50 cytokine/TLR genes
plot_heatmap(
 data_pathways[top_var_pathway_genes, all_samples],
 annotation_col,
 "Top 50 Variable Cytokine (IL, TNF) and Toll-like Receptor Genes"
)
