import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# --------------------------------------------------
# Input files (already extracted)
# --------------------------------------------------
tarbase_file = "Homo_sapiens.tsv"
mirtarbase_file = "hsa_MTI.csv"

hub_genes = ["CSNK2A1", "SUPT16H", "PHIP", "CHD3", "AGO2"]

# --------------------------------------------------
# 1. Load TarBase v9 (PLAIN TSV, NOT gzip)
# --------------------------------------------------
tarbase_df = pd.read_csv(tarbase_file, sep="\t")

# Auto-detect columns
tar_mirna_col = [c for c in tarbase_df.columns if "mirna" in c.lower()][0]
tar_gene_col  = [c for c in tarbase_df.columns if "gene" in c.lower()][0]

# Filter for hub genes
tarbase_hits = tarbase_df[
    tarbase_df[tar_gene_col].isin(hub_genes)
][[tar_gene_col, tar_mirna_col]].dropna().drop_duplicates()

tarbase_hits.columns = ["Gene", "miRNA"]

# --------------------------------------------------
# 2. Load miRTarBase / MTI file
# --------------------------------------------------
mti_df = pd.read_csv(mirtarbase_file)

mti_mirna_col = [c for c in mti_df.columns if "mirna" in c.lower()][0]
mti_gene_col  = [c for c in mti_df.columns if "gene" in c.lower()][0]

mti_hits = mti_df[
    mti_df[mti_gene_col].isin(hub_genes)
][[mti_gene_col, mti_mirna_col]].dropna().drop_duplicates()

mti_hits.columns = ["Gene", "miRNA"]

# --------------------------------------------------
# 3. Save interaction tables
# --------------------------------------------------
tarbase_hits.to_csv(
    "TarBase_hubgene_miRNA_interactions.csv", index=False
)
mti_hits.to_csv(
    "miRTarBase_hubgene_miRNA_interactions.csv", index=False
)

# --------------------------------------------------
# 4. Create miRNA sets
# --------------------------------------------------
tarbase_mirnas = set(tarbase_hits["miRNA"])
mirtarbase_mirnas = set(mti_hits["miRNA"])

shared = tarbase_mirnas & mirtarbase_mirnas
tarbase_only = tarbase_mirnas - mirtarbase_mirnas
mirtarbase_only = mirtarbase_mirnas - tarbase_mirnas

# --------------------------------------------------
# 5. Save miRNA lists
# --------------------------------------------------
pd.Series(sorted(shared), name="miRNA").to_csv(
    "Shared_miRNAs.txt", index=False
)

pd.Series(sorted(tarbase_only), name="miRNA").to_csv(
    "TarBase_only_miRNAs.txt", index=False
)

pd.Series(sorted(mirtarbase_only), name="miRNA").to_csv(
    "miRTarBase_only_miRNAs.txt", index=False
)

# --------------------------------------------------
# 6. Plot Venn diagram
# --------------------------------------------------
plt.figure(figsize=(6, 6))

venn2(
    subsets=(len(tarbase_only), len(mirtarbase_only), len(shared)),
    set_labels=("TarBase v9", "miRTarBase")
)

plt.title(
    "Overlap of miRNA Regulators of Hub Genes\n"
    "(CSNK2A1, SUPT16H, PHIP, CHD3, AGO2)"
)

plt.tight_layout()
plt.savefig(
    "synVenn_TarBase_vs_miRTarBase_HubGene_miRNAs.png",
    dpi=300
)
plt.close()

# --------------------------------------------------
# 7. Summary
# --------------------------------------------------
print("Analysis completed successfully.")
print(f"TarBase-only miRNAs    : {len(tarbase_only)}")
print(f"miRTarBase-only miRNAs : {len(mirtarbase_only)}")
print(f"Shared miRNAs          : {len(shared)}")