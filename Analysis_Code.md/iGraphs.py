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

# Step 1: Define TF dictionaries for each group
il_tf_dict = {
   "IL-9": ["STAT5A", "STAT5B", "IRF4", "BATF", "NFATC2", "RELA"],
   "IL-21": ["STAT3", "BATF", "IRF4", "BCL6", "RELA", "NFATC1"]
}

tlr_tf_dict = {
   "TLR9": ["IRF7", "NFYA", "NFKB1", "SPI1", "IRF5", "RELA", "IRF8", "STAT1", "STAT2", "RUNX1"],
   "TLR7_8": ["IRF5", "IRF7", "NFKB1", "RELA", "STAT1", "IRF4", "FOXP3", "JUN", "FOS", "ELK1"],
   "TLR3": ["IRF3", "IRF7", "NFKB1", "RELA", "STAT2", "IRF1", "ATF2", "FOS", "BATF", "NFATC2"],
   "MyD88_independent_TLR4": ["IRF3", "IRF7", "STAT1", "NFKB1", "RELA"],
   "MyD88_dependent_endosome": ["IRF5", "IRF7", "RELA", "NFKB1", "STAT1"]
}

other_tnf_tf_dict = {
   "Cellular_Responses_to_Stimuli": ["JUN", "FOS", "CREB1", "STAT3", "RELA"],
   "Cellular_Responses_to_Stress": ["FOXO3", "ATF2", "STAT1", "IRF1"],
   "Signal_Transduction": ["NFKB1", "RELA", "STAT3", "ELK1"],
   "Oncogene_Induced_Senescence": ["STAT3", "RELA", "FOXO3", "BCL6"]
}

# Step 2: Load edge file
edges_df = pd.read_csv("top_50%_shared_edges.csv")

# Step 3: Helper functions
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
       print(f"No edges for {title}")
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

# Step 4: Create directories
os.makedirs("il", exist_ok=True)
os.makedirs("tlr", exist_ok=True)
os.makedirs("other_tnf", exist_ok=True)

# Step 5: Plot all networks
for cat_name, tf_dict, folder in [
   ("IL", il_tf_dict, "il"),
   ("TLR", tlr_tf_dict, "tlr"),
   ("Other TNF", other_tnf_tf_dict, "other_tnf")
]:
   for pathway, tf_list in tf_dict.items():
       sub_edges = edges_df[
           (edges_df["X.START_ID"].isin(tf_list)) & (edges_df["X.END_ID"].isin(tf_list))
       ]
       safe_name = f"{folder}/{pathway.replace('/', '_')}_network.png"
       plot_tf_network_weighted_pretty(
           sub_edges,
           tf_list,
           f"{pathway} Network",
           safe_name
       )
