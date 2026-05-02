import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import math
import numpy as np

# --------------------------------------------------
# Input
# --------------------------------------------------
edges_file = "BrainSpecific_miRNA_Gene_Interactions_RPM15.csv"

hub_genes = [
    "CSNK2A1",
    "SUPT16H",
    "PHIP",
    "CHD3",
    "AGO2"
]

master_mirnas = [
    "hsa-let-7b-5p",
    "hsa-miR-186-5p",
    "hsa-miR-16-5p",
    "hsa-let-7a-5p",
    "hsa-miR-15a-5p",
    "hsa-miR-299-5p",
    "hsa-miR-218-5p",
    "hsa-miR-21-5p",
    "hsa-miR-195-5p",
    "hsa-miR-15b-5p",
    "hsa-miR-320b",
    "hsa-let-7c-5p",
    "hsa-miR-103a-3p",
    "hsa-miR-107"
]

# --------------------------------------------------
# Load & filter data
# --------------------------------------------------
edges = pd.read_csv(edges_file)
edges = edges[edges["Gene"].isin(hub_genes)]

genes = sorted(edges["Gene"].unique())
mirnas = sorted(edges["miRNA"].unique())

master_nodes = [m for m in mirnas if m in master_mirnas]
other_mirnas = [m for m in mirnas if m not in master_nodes]

# --------------------------------------------------
# Graph
# --------------------------------------------------
G = nx.Graph()
for _, r in edges.iterrows():
    G.add_edge(r["miRNA"], r["Gene"])

# --------------------------------------------------
# MULTI-RING LAYOUT (NO OVERLAP)
# --------------------------------------------------
pos = {}

# Hub genes (inner ring)
gene_radius = 2.0
for i, g in enumerate(genes):
    angle = 2 * math.pi * i / len(genes)
    pos[g] = (gene_radius * math.cos(angle),
              gene_radius * math.sin(angle))

# Master miRNAs (ring 2)
master_radius = 4.0
for i, m in enumerate(master_nodes):
    angle = 2 * math.pi * i / len(master_nodes)
    pos[m] = (master_radius * math.cos(angle),
              master_radius * math.sin(angle))

# Other miRNAs → split into TWO rings
ring1 = other_mirnas[:len(other_mirnas)//2]
ring2 = other_mirnas[len(other_mirnas)//2:]

outer_radius_1 = 6.5
outer_radius_2 = 8.5

for i, m in enumerate(ring1):
    angle = 2 * math.pi * i / len(ring1)
    pos[m] = (outer_radius_1 * math.cos(angle),
              outer_radius_1 * math.sin(angle))

for i, m in enumerate(ring2):
    angle = 2 * math.pi * i / len(ring2)
    pos[m] = (outer_radius_2 * math.cos(angle),
              outer_radius_2 * math.sin(angle))

# --------------------------------------------------
# Plot
# --------------------------------------------------
plt.figure(figsize=(30, 30))

nx.draw_networkx_nodes(G, pos, nodelist=genes,
                       node_color="#f08080",
                       node_size=5600,
                       edgecolors="black",
                       linewidths=3,
                       label="Hub genes")

nx.draw_networkx_nodes(G, pos, nodelist=master_nodes,
                       node_color="#0b3cde",
                       node_size=3400,
                       edgecolors="black",
                       linewidths=2.5,
                       label="Master regulator miRNAs")

nx.draw_networkx_nodes(G, pos, nodelist=other_mirnas,
                       node_color="#87cefa",
                       node_size=2600,
                       edgecolors="black",
                       linewidths=1.8,
                       label="Other miRNAs")

nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.45)

# --------------------------------------------------
# Labels (offset by node type)
# --------------------------------------------------
label_pos = {}
for n in pos:
    if n in genes:
        label_pos[n] = (pos[n][0], pos[n][1] - 0.6)
    elif n in master_nodes:
        label_pos[n] = (pos[n][0], pos[n][1] - 0.4)
    else:
        label_pos[n] = (pos[n][0], pos[n][1] - 0.35)

nx.draw_networkx_labels(G, label_pos,
                        font_size=15,
                        font_weight="bold")

# --------------------------------------------------
# Title (unchanged)
# --------------------------------------------------
plt.text(
    0.5, -0.05,
    "Brain-Specific miRNA–Hub Gene Regulatory Network of Syndromic Autism Hubgenes (RPM ≥ 10)\n"
    "Master regulator miRNAs highlighted in dark blue",
    fontsize=24,
    fontweight="bold",
    ha="center",
    transform=plt.gca().transAxes
)

plt.legend(loc="upper right", bbox_to_anchor=(0.98, 0.98),
           fontsize=18, frameon=True)

plt.axis("off")
plt.tight_layout()

plt.savefig(
    "BrainSpecific_miRNA_HubGene_Network_SYNDROMIC_FINAL.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Syndromic autism network generated with multi-ring layout (no overlaps).")
