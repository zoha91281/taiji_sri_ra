!pip install igraph cairocffi

from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from igraph import Graph, plot
import numpy as np
import pandas as pd
from IPython.display import Image, display
from google.colab import files
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# Load PageRank scores
pagerank_df = pd.read_csv("PageRanks_updated_names.csv")
pagerank_df.columns = [c.strip() for c in pagerank_df.columns]

# Compute average PageRank per cohort for each TF
pagerank_df['avg_OA'] = pagerank_df.loc[:, [c for c in pagerank_df.columns if c.startswith('OA')]].mean(axis=1)
pagerank_df['avg_RA'] = pagerank_df.loc[:, [c for c in pagerank_df.columns if c.startswith('RA')]].mean(axis=1)

# Create dicts for quick lookup: TF -> avg PageRank per cohort
pagerank_OA_dict = dict(zip(pagerank_df['V1'], pagerank_df['avg_OA']))
pagerank_RA_dict = dict(zip(pagerank_df['V1'], pagerank_df['avg_RA']))

# TLR TF dictionary (only TLR9 and TLR7_8)
tlr_tf_dict = {
    "TLR9": ["IRF7", "NFYA", "NFKB1", "SPI1", "IRF5", "RELA", "IRF8", "STAT1", "STAT2", "RUNX1"],
    "TLR7_8": ["IRF5", "IRF7", "NFKB1", "RELA", "STAT1", "IRF4", "FOXP3", "JUN", "FOS", "ELK1"]
}

# Load edges CSV
edges_df = pd.read_csv("top_50%_shared_edges (1).csv")  # update filename if needed
print("Unique sample labels:", edges_df['sample'].unique())

# Create output folder
os.makedirs("tlr", exist_ok=True)

def add_isolated_nodes(g, all_nodes):
    present_nodes = set(g.vs['label'])
    missing_nodes = set(all_nodes) - present_nodes
    for node in missing_nodes:
        g.add_vertex(name=node, label=node, size=15, color="#4a90e2")
    return g

def plot_tf_network_weighted_pretty(df_edges, all_tfs, title, filename,
                                   pagerank_dict_cohort,
                                   color_start="#74add1", color_end="#d73027",
                                   compact_factor=0.7):
    if df_edges.empty:
        print(f"⚠️ No edges for {title}")
        return

    unique_nodes = pd.unique(df_edges[['X.START_ID', 'X.END_ID']].values.ravel())
    node_index = {node: i for i, node in enumerate(unique_nodes)}
    edge_list = [(node_index[s], node_index[t]) for s, t in df_edges[['X.START_ID', 'X.END_ID']].values]
    weights = df_edges['weight'].tolist()

    g = Graph(edge_list, directed=False)
    g.simplify(multiple=True, loops=True)
    g.vs['label'] = unique_nodes
    g.es['weight'] = weights

    g = add_isolated_nodes(g, all_tfs)

    # Use cohort-specific PageRank for node size and color
    pagerank_scores = np.array([pagerank_dict_cohort.get(v["label"], 0.01) for v in g.vs])
    pr_min, pr_max = pagerank_scores.min(), pagerank_scores.max()
    norm_pageranks = (pagerank_scores - pr_min) / (pr_max - pr_min + 1e-9)

    node_sizes = 30 + norm_pageranks * 50
    g.vs['size'] = node_sizes.tolist()

    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", [color_start, color_end])
    node_colors_rgb = cmap(norm_pageranks)
    node_colors_hex = [mcolors.to_hex(c) for c in node_colors_rgb]

    for v in g.vs:
        if v['label'] not in unique_nodes:
            v['color'] = "#4a90e2"
        else:
            v['color'] = node_colors_hex[v.index]

    layout = g.layout_kamada_kawai()
    coords = np.array(layout.coords)
    coords *= compact_factor
    layout = coords.tolist()

    scaled_weights = np.interp(weights, (min(weights), max(weights)), (1.0, 3.0))
    edge_color = (0.6, 0.6, 0.6, 0.4)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_facecolor("#fafafa")
    for spine in ax.spines.values():
        spine.set_edgecolor('#dcdcdc')
        spine.set_linewidth(1.0)

    ax.set_title(title, fontsize=14, fontweight='normal', pad=10)

    plot(
        g,
        target=ax,
        layout=layout,
        vertex_label=g.vs['label'],
        vertex_color=g.vs['color'],
        vertex_size=g.vs['size'],
        vertex_frame_color=None,
        vertex_label_size=10,
        vertex_label_dist=0.5,
        vertex_label_color='black',
        edge_width=scaled_weights,
        edge_color=edge_color,
        margin=30,
        bbox=(600, 600),
        background="#fafafa"
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="4%", pad=0.3)
    cb = ColorbarBase(cax, cmap=cmap, norm=mcolors.Normalize(vmin=pr_min, vmax=pr_max), orientation='horizontal')
    cb.set_label('')
    cax.xaxis.set_ticks_position('bottom')
    cax.xaxis.set_label_position('bottom')
    cax.text(0.5, 1.25, 'Node PageRank (regulatory importance)', ha='center', va='bottom', fontsize=10, transform=cax.transAxes)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(filename, dpi=300)
    plt.close(fig)

    display(Image(filename=filename))
    files.download(filename)
    print(f"✅ {title} rendered and saved as {filename}")

# Settings
top_n = 25     # Top edges
top_m_tfs = 8  # Top TFs per pathway by PageRank

selected_pathways = ["TLR9", "TLR7_8"]

for pathway in selected_pathways:
    tf_list = tlr_tf_dict[pathway]
    print(f"\n🔬 {pathway} pathway")

    for cohort in ["OA", "RA"]:
        print(f"Processing cohort: {cohort}")

        # Choose appropriate PageRank dict by cohort
        pagerank_dict_cohort = pagerank_OA_dict if cohort == "OA" else pagerank_RA_dict

        # Rank TFs by cohort-specific PageRank
        tf_ranks = [(tf, pagerank_dict_cohort.get(tf, 0)) for tf in tf_list]
        tf_ranks = sorted(tf_ranks, key=lambda x: x[1], reverse=True)
        top_tfs = [tf for tf, _ in tf_ranks[:top_m_tfs]]

        cohort_edges = edges_df[edges_df['sample'].str.lower().str.startswith(cohort.lower())].copy()

        # Filter edges to those between top TFs only
        edges_touching_tfs = cohort_edges[
            (cohort_edges["X.START_ID"].isin(top_tfs)) |
            (cohort_edges["X.END_ID"].isin(top_tfs))
        ].copy()

        sub_edges = edges_touching_tfs[
            (edges_touching_tfs["X.START_ID"].isin(top_tfs)) &
            (edges_touching_tfs["X.END_ID"].isin(top_tfs))
        ].copy()

        sub_edges = sub_edges.sort_values('weight', ascending=False).head(top_n)

        # Symmetrize for undirected graph
        reversed_edges = sub_edges.rename(columns={"X.START_ID": "X.END_ID", "X.END_ID": "X.START_ID"})
        sym_edges = pd.concat([sub_edges, reversed_edges]).drop_duplicates()

        title = f"{pathway} Network ({cohort})"
        filename = f"tlr/{pathway}_{cohort}.png"

        plot_tf_network_weighted_pretty(sym_edges, top_tfs, title, filename, pagerank_dict_cohort)
