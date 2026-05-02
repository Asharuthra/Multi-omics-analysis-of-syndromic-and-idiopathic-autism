#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================
# User Configuration
# ==============================
HPA_DIR = "/home/nivi/Downloads/autism_brainspan/tsv-file"  # folder with per-gene TSVs
OUTPUT_DIR = "outputs_hpa"
GENES_OF_INTEREST = ["UBA52", "RPS27A", "KAT2B", "AR", "H3C1"]

# Match GTEx-style brain order (excluding spinal cord)
BRAIN_REGIONS = [
    "amygdala",
    "basal ganglia",
    "cerebellum",
    "cerebral cortex",
    "choroid plexus",
    "hippocampal formation",
    "hypothalamus",
    "medulla oblongata",
    "midbrain",
    "pons",
    "thalamus",
    "white matter"
]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# Load all per-gene TSVs (HPA data)
# ==============================
all_data = {}
for fname in os.listdir(HPA_DIR):
    if not fname.endswith(".tsv"):
        continue
    path = os.path.join(HPA_DIR, fname)
    gene = os.path.splitext(fname)[0].strip()
    if gene not in GENES_OF_INTEREST:
        continue

    try:
        df = pd.read_csv(path, sep="\t")
        df.columns = df.columns.str.strip()
        expr_cols = [c for c in df.columns if "Brain RNA" in c and "nTPM" in c]
        if not expr_cols:
            print(f"⚠️ No nTPM columns found for {gene}")
            continue
        values = df.loc[0, expr_cols].astype(float)
        all_data[gene] = values
    except Exception as e:
        print(f"⚠️ Error reading {fname}: {e}")

if not all_data:
    raise ValueError("❌ No valid gene expression data found in TSV folder!")

# ==============================
# Combine into DataFrame
# ==============================
expr = pd.DataFrame(all_data).T
expr.columns = [c.replace("Brain RNA - ", "").replace(" [nTPM]", "") for c in expr.columns]

# Normalize brain region column order (exclude spinal cord)
expr = expr[[col for col in BRAIN_REGIONS if col in expr.columns]]

# ==============================
# Log2 Transform
# ==============================
df = np.log2(expr + 1).replace([np.inf, -np.inf], 0).fillna(0)

# ==============================
# Save Processed CSV
# ==============================
csv_path = os.path.join(OUTPUT_DIR, "hpa_brain_expression.csv")
df.to_csv(csv_path)

# ==============================
# Publication-Ready Heatmap
# ==============================
plt.figure(figsize=(14, 4))
sns.heatmap(
    df,
    cmap="cividis",
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "log2(nTPM + 1)"}
)

plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=12, weight="bold")

plt.title("Brain Regional Expression Patterns of Idiopathic Autism Hub Genes (via HPA)", fontsize=16, weight="bold")	
plt.tight_layout()

heatmap_path = os.path.join(OUTPUT_DIR, "hpa_heatmap_publication.png")
plt.savefig(heatmap_path, dpi=600)
plt.close()

print(f"✅ Saved CSV: {csv_path}")
print(f"✅ Saved heatmap: {heatmap_path}")

