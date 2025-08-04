import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import gseapy as gp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

# Step 1: Load your data and preprocess
df = pd.read_csv("PageRanks_updated_names.csv", sep=",")

oa_cols = [col for col in df.columns if col.startswith("OA")]
ra_cols = [col for col in df.columns if col.startswith("RA")]

df[oa_cols] = df[oa_cols].apply(pd.to_numeric, errors='coerce')
df[ra_cols] = df[ra_cols].apply(pd.to_numeric, errors='coerce')

df['mean_OA'] = df[oa_cols].mean(axis=1)
df['mean_RA'] = df[ra_cols].mean(axis=1)
df['log2FC'] = np.log2(df['mean_RA'] + 1e-9) - np.log2(df['mean_OA'] + 1e-9)

# Calculate p-values per TF
pvals = []
for _, row in df.iterrows():
   ra_vals = row[ra_cols].dropna().values.astype(float)
   oa_vals = row[oa_cols].dropna().values.astype(float)
   if len(ra_vals) > 1 and len(oa_vals) > 1:
       _, p = ttest_ind(ra_vals, oa_vals, equal_var=False)
   else:
       p = 1.0
   pvals.append(p)
df['pval'] = pvals

fc_threshold = np.log2(1.2)
filtered_tfs = df[(df['pval'] < 0.05) & (np.abs(df['log2FC']) > fc_threshold)]

print(f"Total TFs passing filters: {len(filtered_tfs)}")
print(f"Upregulated TFs: {(filtered_tfs['log2FC'] > fc_threshold).sum()}")
print(f"Downregulated TFs: {(filtered_tfs['log2FC'] < -fc_threshold).sum()}")

# Step 2: Run Reactome pathway enrichment using Enrichr API
tf_list = filtered_tfs['V1'].dropna().tolist()

enr = gp.enrichr(gene_list=tf_list,
                gene_sets='Reactome_2022',
                organism='Human',
                outdir=None,
                cutoff=1.0)

reactome_res = enr.results

# Filter pathways with adjusted p-value (FDR) < 0.05
sig_pathways = reactome_res[reactome_res['Adjusted P-value'] < 0.01]

print(f"Significant Reactome pathways found: {len(sig_pathways)}")

# Step 3: Prepare data for dotplot
# Use Combined Score as "NES"-like value (scaled influence)
nes_values = sig_pathways['Combined Score'].to_numpy(dtype=float)
terms = sig_pathways['Term']
pvals = sig_pathways['Adjusted P-value'].to_numpy(dtype=float)

dot_size = -np.log10(pvals) * 100  # size scaled by significance
norm = Normalize(vmin=nes_values.min(), vmax=nes_values.max())

custom_cmap = LinearSegmentedColormap.from_list("blue_to_red", ["#74add1", "#d73027"])
colors = custom_cmap(norm(nes_values))

# Step 4: Plot dotplot
fig, ax = plt.subplots(figsize=(10, 15))
scatter = ax.scatter(nes_values, terms, s=dot_size, c=colors)

ax.set_xlabel('Combined Score (Reactome Pathway Enrichment)')
ax.set_ylabel('Reactome Pathways')
ax.set_title('Reactome Pathway Enrichment Dotplot (Filtered TFs)')

sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Combined Score')

plt.tight_layout()
# Save plot as PNG
plt.savefig("reactome_enrichment_dotplot.png", dpi=300, bbox_inches='tight')

# Optionally display plot in interactive session
plt.show()
