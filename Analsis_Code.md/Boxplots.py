import pandas as pd
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt

# Parameters
pagerank_path = "rnaseq_updated_names.csv"  # Your expression matrix
wilcox_path = "wilcoxon_results_final (1).csv"  # Wilcoxon DE results
log2fc_thresh = 0.5

# -------- Load and prepare pagerank/expression data --------
pagerank_df = pd.read_csv(pagerank_path)

# Use 'gene_name' as index, drop 'gene_id' column (which has NaNs)
pagerank_df = pagerank_df.set_index('gene_name').drop(columns=['gene_id'])

# Transpose so samples are rows, genes are columns
pagerank_df = pagerank_df.T

# Convert all values to numeric (force errors to NaN)
pagerank_df = pagerank_df.apply(pd.to_numeric, errors='coerce')

print("Pagerank shape after cleanup:", pagerank_df.shape)
print(pagerank_df.head())
print(pagerank_df.dtypes)

# -------- Z-score scaling per gene --------
scaler = StandardScaler()
pagerank_scaled_values = scaler.fit_transform(pagerank_df)
pagerank_scaled = pd.DataFrame(pagerank_scaled_values,
                               columns=pagerank_df.columns,
                               index=pagerank_df.index)

# Add group label (samples starting with 'OA' or 'RA')
pagerank_scaled["Group"] = pagerank_scaled.index.to_series().apply(lambda x: "OA" if x.startswith("OA") else "RA")

# -------- Load and prepare Wilcoxon DE results --------
wilcox_df = pd.read_csv(wilcox_path)

# Rename log2FC column correctly
if 'log2FC_RA_vs_OA' in wilcox_df.columns:
    wilcox_df.rename(columns={'log2FC_RA_vs_OA': 'log2FC'}, inplace=True)

# Create abs_log2FC if missing
if 'abs_log2FC' not in wilcox_df.columns and 'log2FC' in wilcox_df.columns:
    wilcox_df['abs_log2FC'] = wilcox_df['log2FC'].abs()

# Uppercase gene names for matching
wilcox_df['gene_name'] = wilcox_df['gene_name'].astype(str).str.upper()

# -------- Function to get top 5 TFs by pathway prefix --------
def get_top5_tfs_by_pathway(prefix):
    subset = wilcox_df[wilcox_df['gene_name'].str.startswith(prefix)]
    filtered = subset[(subset['p_value'] < 0.05) & (subset['abs_log2FC'] > log2fc_thresh)]
    if len(filtered) < 5:
        filtered = subset[subset['abs_log2FC'] > log2fc_thresh]
    if len(filtered) < 5:
        filtered = subset.sort_values(by=['p_value', 'abs_log2FC'], ascending=[True, False]).head(5)
    else:
        filtered = filtered.sort_values(by=['p_value', 'abs_log2FC'], ascending=[True, False]).head(5)
    return filtered['gene_name'].tolist()

# -------- Get top 5 genes for each category --------
top_il = get_top5_tfs_by_pathway('IL')
top_tnf = get_top5_tfs_by_pathway('TNF')
top_tlr = get_top5_tfs_by_pathway('TLR')

# Top 5 overall DEGs (filtered and sorted)
all_degs = wilcox_df[(wilcox_df['p_value'] < 0.05) & (wilcox_df['abs_log2FC'] > log2fc_thresh)]
top_all_degs = all_degs.sort_values(by=['p_value', 'abs_log2FC'], ascending=[True, False]).head(5)['gene_name'].tolist()

print(f"Top 5 all DEGs: {top_all_degs}")
print(f"Top 5 IL: {top_il}")
print(f"Top 5 TNF: {top_tnf}")
print(f"Top 5 TLR: {top_tlr}")

# -------- Prepare data for plotting --------
pagerank_long = pagerank_scaled.reset_index().melt(
    id_vars=["index", "Group"], var_name="Gene", value_name="Zscore"
)
pagerank_long.rename(columns={"index": "Sample"}, inplace=True)

# -------- Plotting function --------
def plot_boxplot(genes, title, filename):
    subset = pagerank_long[pagerank_long["Gene"].isin(genes)]
    if subset.empty:
        print(f"WARNING: No data to plot for genes: {genes}")
        return
    plt.figure(figsize=(12, 7))
    sns.boxplot(data=subset, x="Gene", y="Zscore", hue="Group",
                palette={"OA": "#74add1", "RA": "#d73027"})
    plt.title(title, fontsize=16)
    plt.xlabel("Gene/TF", fontsize=14)
    plt.ylabel("Z-scored Expression/PageRank", fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()

# -------- Plot all groups --------
plot_boxplot(top_all_degs, "Top 5 All DEGs", "boxplot_all_degs.png")
plot_boxplot(top_il, "Top 5 IL Genes", "boxplot_il.png")
plot_boxplot(top_tnf, "Top 5 TNF Genes", "boxplot_tnf.png")
plot_boxplot(top_tlr, "Top 5 TLR Genes", "boxplot_tlr.png")
