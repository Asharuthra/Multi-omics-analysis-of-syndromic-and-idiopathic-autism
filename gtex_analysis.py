#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================
# User Configuration
# ==============================
GCT_PATH = "/home/nivi/Downloads/autism_brainspan/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct"
META_PATH = "/home/nivi/Downloads/autism_brainspan/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
OUTPUT_DIR = "outputs"
GENES_OF_INTEREST = ["UBA52", "RPS27A", "KAT2B", "AR", "H3C1"]

BRAIN_REGIONS = [
    "Brain - Amygdala",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Caudate (basal ganglia)",
    "Brain - Cerebellar Hemisphere",
    "Brain - Cerebellum",
    "Brain - Cortex",
    "Brain - Frontal Cortex (BA9)",
    "Brain - Hippocampus",
    "Brain - Hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Substantia nigra"
]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# Load Metadata (Brain Regions Only)
# ==============================
meta = pd.read_csv(META_PATH, sep="\t", dtype=str)
meta = meta[(meta["SMTS"] == "Brain") & (meta["SMTSD"].isin(BRAIN_REGIONS))]
sample_to_region = meta.set_index("SAMPID")["SMTSD"]

# ==============================
# STREAM-LOAD GCT FILE (NO FULL READ)
# ==============================
needed_samples = set(sample_to_region.index)
genes_needed = set(GENES_OF_INTEREST)
gene_data = {}

with open(GCT_PATH, "r") as f:
    next(f); next(f)  # Skip 2 metadata lines
    header = next(f).strip().split("\t")

    col_indices = [i for i, col in enumerate(header) if col in needed_samples]
    sample_names = [header[i] for i in col_indices]

    for line in f:
        parts = line.strip().split("\t")
        gene_symbol = parts[1]  # Description column

        if gene_symbol in genes_needed:
            values = [float(parts[i]) for i in col_indices]
            gene_data[gene_symbol] = dict(zip(sample_names, values))

expr = pd.DataFrame(gene_data).T  # Genes x Samples

# ==============================
# Aggregate by Brain Region
# ==============================
expr = expr.T
expr["Region"] = expr.index.map(sample_to_region)
df = expr.groupby("Region").mean().T  # Genes x Regions

# ==============================
# Log2 Transform
# ==============================
df = np.log2(df + 1).replace([np.inf, -np.inf], 0).fillna(0)
df = df[BRAIN_REGIONS]  # Keep ordered regions

# ==============================
# Publication-Ready Heatmap (PNG)
# ==============================
plt.figure(figsize=(14, 4))
sns.heatmap(
    df,
    cmap="cividis",
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "log2(TPM + 1)"}
)

plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=12, weight='bold')

plt.title("GTEx Brain Expression of Idiopathic Autism Linked Hub Genes", fontsize=16, weight='bold')
plt.tight_layout()

heatmap_path = os.path.join(OUTPUT_DIR, "gtex_heatmap_publication.png")
plt.savefig(heatmap_path, dpi=600)
plt.close()

csv_path = os.path.join(OUTPUT_DIR, "gtex_brain_expression.csv")
df.to_csv(csv_path)

print(f"✅ Saved CSV: {csv_path}")
print(f"✅ Saved heatmap: {heatmap_path}")

