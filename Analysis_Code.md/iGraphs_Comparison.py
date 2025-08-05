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

# Define TLR TF dictionaries only
tlr_tf_dict = {
   "TLR9": ["IRF7", "NFYA", "NFKB1", "SPI1", "IRF5", "RELA", "IRF8", "STAT1", "STAT2", "RUNX1"],
   "TLR7_8": ["IRF5", "IRF7", "NFKB1", "RELA", "STAT1", "IRF4", "FOXP3", "JUN", "FOS", "ELK1"],
   "TLR3": ["IRF3", "IRF7", "NFKB1", "RELA", "STAT2", "IRF1", "ATF2", "FOS", "BATF", "NFATC2"],
   "MyD88_independent_TLR4": ["IRF3", "IRF7", "STAT1", "NFKB1", "RELA"],
   "MyD88_dependent_endosome": ["IRF5", "IRF7", "RELA", "NFKB1", "STAT1"]
}

# Load your edges CSV (update filename as needed)
edges_df = pd.read_csv("top_50%_shared_edges (1).csv")

print("Unique sample labels:", edges_df['sample'].unique())

def add_isolated_nodes(g, all_nodes):
    present_nodes = set(g.vs['label'])
    missing_nodes = set(all_nodes) - present_nodes
    for node in missing_nodes:
        g.add_vertex(name=node, label=node, size=15, color="#4a90e2")
    return g

def plot_tf_network_weighted_pretty(df_edges, all_tfs, title, filename,
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

    degrees = np.array(g.degree())
    deg_min, deg_max = degrees.min(), degrees.max()
    norm_degrees = (degrees - deg_min) / (deg_max - deg_min + 1e-9)
    node_sizes = 30 + norm_degrees * 50
    g.vs['size'] = node_sizes.tolist()

    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", [color_start, color_end])
    node_colors_rgb = cmap(norm_degrees)
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
    cb = ColorbarBase(cax, cmap=cmap, norm=mcolors.Normalize(vmin=deg_min, vmax=deg_max), orientation='horizontal')
    cb.set_label('')
    cax.xaxis.set_ticks_position('bottom')
    cax.xaxis.set_label_position('bottom')
    cax.text(0.5, 1.25, 'Node degree (number of connections)', ha='center', va='bottom', fontsize=10, transform=cax.transAxes)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(filename, dpi=300)
    plt.close(fig)

    display(Image(filename=filename))
    files.download(filename)
    print(f"✅ {title} rendered and saved as {filename}")

# Output folder
os.makedirs("tlr", exist_ok=True)

def get_cohort(sample):
    return "OA" if sample.lower().startswith("oa") else "RA"

top_n = 25  # Change to 10 if you want only top 10 edges

for pathway, tf_list in tlr_tf_dict.items():
    print(f"\n🔬 {pathway} pathways")
    for cohort in ["OA", "RA"]:
        print(f"Processing cohort: {cohort}")
        
        cohort_edges = edges_df[edges_df['sample'].str.lower().str.startswith(cohort.lower())].copy()
        
        # Filter edges touching any TF in this pathway (union)
        edges_touching_tfs = cohort_edges[
            (cohort_edges["X.START_ID"].isin(tf_list)) | (cohort_edges["X.END_ID"].isin(tf_list))
        ].copy()
        
        # Keep only edges between TFs in pathway (intersection)
        sub_edges = edges_touching_tfs[
            (edges_touching_tfs["X.START_ID"].isin(tf_list)) &
            (edges_touching_tfs["X.END_ID"].isin(tf_list))
        ].copy()

        # Limit to top_n edges by weight
        sub_edges = sub_edges.sort_values('weight', ascending=False).head(top_n)

        # Symmetrize edges to ensure undirected graph completeness
        reversed_edges = sub_edges.rename(columns={"X.START_ID": "X.END_ID", "X.END_ID": "X.START_ID"})
        sym_edges = pd.concat([sub_edges, reversed_edges]).drop_duplicates()

        title = f"{pathway} Network ({cohort})"
        filename = f"tlr/{pathway}_{cohort}.png"

        plot_tf_network_weighted_pretty(sym_edges, tf_list, title, filename)
