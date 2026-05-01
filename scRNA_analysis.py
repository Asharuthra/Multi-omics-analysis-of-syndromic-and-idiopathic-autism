#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import anndata
import matplotlib.pyplot as plt
import seaborn as sns

# ======================================
# USER CONFIGURATION
# ======================================
DATA_DIR = "/home/nivi/Downloads/autism_brainspan/allen"
EXPR_FILE = os.path.join(DATA_DIR, "human_MTG_2018-06-14_exon-matrix.csv")
GENE_FILE = os.path.join(DATA_DIR, "human_MTG_2018-06-14_genes-rows.csv")
META_FILE = os.path.join(DATA_DIR, "human_MTG_2018-06-14_samples-columns.csv")
OUT_DIR = os.path.join(DATA_DIR, "outputs_final")
os.makedirs(OUT_DIR, exist_ok=True)

# syndromic autism hub genes
GENES = ["CSNK2A1", "SUPT16H", "PHIP", "CHD3", "AGO2"]

# ======================================
# STEP 1. LOAD METADATA
# ======================================
print("📂 Loading metadata...")
meta = pd.read_csv(META_FILE)
meta.index = meta.iloc[:, 0].astype(str)

genes = pd.read_csv(GENE_FILE)
gene_col_candidates = ["gene_symbol", "gene", "gene_name", "external_gene_name", "symbol"]
gene_col = next((c for c in gene_col_candidates if c in genes.columns), genes.columns[0])
genes.index = genes[gene_col].astype(str)
print(f"✅ Using gene column: {gene_col}")

# ======================================
# STEP 2. STREAM EXPRESSION MATRIX (LOW MEMORY)
# ======================================
print("🔄 Streaming expression CSV in chunks...")
reader = pd.read_csv(EXPR_FILE, chunksize=500, index_col=0)
chunks = []
row_count = 0

for chunk in reader:
    chunk = chunk.apply(pd.to_numeric, errors='coerce').fillna(0)
    chunk_sparse = sparse.csr_matrix(chunk.values.astype(np.float32))
    chunks.append(chunk_sparse)
    row_count += chunk.shape[0]
    if row_count % 5000 == 0:
        print(f"  processed {row_count} genes...")

expr_sparse = sparse.vstack(chunks, format='csr')
print(f"✅ Sparse matrix shape: {expr_sparse.shape}")

# ======================================
# STEP 3. TRANSPOSE MATRIX
# ======================================
expr_sparse_T = expr_sparse.T  # cells × genes
expr_genes = genes.iloc[:expr_sparse.shape[0]]
expr_meta = meta.iloc[:expr_sparse_T.shape[0]]
print(f"✅ After transpose: {expr_sparse_T.shape} (cells × genes)")

# ======================================
# STEP 4. BUILD ANNDATA
# ======================================
adata = anndata.AnnData(X=expr_sparse_T, obs=expr_meta, var=expr_genes)
adata.obs_names_make_unique()
adata.var_names_make_unique()

adata_path = os.path.join(OUT_DIR, "human_MTG_sparse_final.h5ad")
adata.write(adata_path, compression="gzip")
print(f"✅ AnnData saved: {adata_path}")

# ======================================
# STEP 5. DEFINE BROAD CELL TYPE MAPPING
# ======================================
print("🔍 Assigning cell type categories...")
broad_map = {
    "Astro": "Astrocyte",
    "Oligo": "Oligodendrocyte",
    "OPC": "OPC",
    "Inh": "Inhibitory neuron",
    "Exc": "Excitatory neuron",
    "Endo": "Endothelial",
    "Micro": "Microglia"
}

adata.obs["broad_type"] = adata.obs["cluster"].apply(
    lambda x: next((broad_map[k] for k in broad_map if k in str(x)), "Other")
)

# Exclude “Other”
adata = adata[adata.obs["broad_type"] != "Other"].copy()
print(f"✅ Remaining cell types: {adata.obs['broad_type'].unique().tolist()}")

# ======================================
# STEP 6. SUBSET HUB GENES AND COMPUTE MEAN EXPRESSION
# ======================================
adata.var[gene_col] = adata.var[gene_col].astype(str)
gene_mask = adata.var[gene_col].isin(GENES)
subset = adata[:, gene_mask]

expr = subset.to_df()
expr["broad_type"] = adata.obs["broad_type"].values
mean_expr = expr.groupby("broad_type").mean().T
mean_expr = np.log2(mean_expr + 1)

# ======================================
# STEP 7. ANOVA + TUKEY POST-HOC
# ======================================
print("📊 Performing one-way ANOVA + Tukey’s HSD...")
anova_records = []
tukey_records = []
significance = {}

for g in GENES:
    data = subset.to_df()[g]
    groups = adata.obs["broad_type"]

    # Perform ANOVA
    grouped_values = [data[groups == ct] for ct in groups.unique()]
    if all(len(v) > 1 for v in grouped_values):
        fval, pval = f_oneway(*grouped_values)
        anova_records.append({"Gene": g, "F_value": fval, "p_value": pval})

        if pval < 0.05:
            df_tukey = pd.DataFrame({"value": data.values, "group": groups.values})
            tukey = pairwise_tukeyhsd(df_tukey["value"], df_tukey["group"], alpha=0.05)
            sig_cells = set()
            for res in tukey.summary().data[1:]:
                g1, g2, meandiff, p_adj, lower, upper, reject = res
                tukey_records.append({
                    "Gene": g,
                    "Group1": g1,
                    "Group2": g2,
                    "MeanDiff": meandiff,
                    "p_adj": p_adj,
                    "Lower": lower,
                    "Upper": upper,
                    "Reject": reject
                })
                if reject:
                    sig_cells.update([g1, g2])
            significance[g] = sig_cells

# Save statistical results
anova_df = pd.DataFrame(anova_records)
anova_path = os.path.join(OUT_DIR, "anova_results_idiopathic.csv")
anova_df.to_csv(anova_path, index=False)

tukey_df = pd.DataFrame(tukey_records)
tukey_path = os.path.join(OUT_DIR, "tukey_results_idiopathic.csv")
tukey_df.to_csv(tukey_path, index=False)

print(f"✅ ANOVA results saved: {anova_path}")
print(f"✅ Tukey post-hoc results saved: {tukey_path}")

# ======================================
# STEP 8. GENERATE ANNOTATED HEATMAP MATRIX
# ======================================
annot_df = pd.DataFrame("", index=mean_expr.index, columns=mean_expr.columns)
for g in mean_expr.index:
    for c in mean_expr.columns:
        val = f"{mean_expr.loc[g, c]:.2f}"
        if (c in significance.get(g, [])) and (mean_expr.loc[g, c] > 1.0):
            val += "*"
        annot_df.loc[g, c] = val

# ======================================
# STEP 9. VISUALIZATION
# ======================================
plt.figure(figsize=(14, 6))
sns.heatmap(
    mean_expr,
    cmap="YlGnBu",
    annot=annot_df,
    fmt="",
    linewidths=0.5,
    cbar_kws={"label": "log2(Expression + 1)"}
)

plt.title(
    "Cell-type Specific Expression of Syndromic Autism Hub Genes\n"
    "(Allen Brain Atlas - Human Middle Temporal Gyrus)",
    fontsize=16,
    weight="bold",
)
plt.xlabel("Broad Cell Type", fontsize=12)
plt.ylabel("Hub Genes", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

out_png = os.path.join(OUT_DIR, "allen_broad_celltype_heatmap_syndromic.png")
plt.savefig(out_png, dpi=600)
plt.close()
print(f"✅ Saved heatmap: {out_png}")

print("🎯 Done — heatmap shows biologically (log2>1) and statistically (*) significant cell-type enrichment.")

