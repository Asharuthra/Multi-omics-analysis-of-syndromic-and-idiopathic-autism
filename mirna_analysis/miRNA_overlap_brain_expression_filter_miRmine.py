import pandas as pd

# --------------------------------------------------
# Input files
# --------------------------------------------------
shared_mirna_file = "Shared_miRNAs.txt"
mirmine_file = "miRmine - Human miRNA Expression Database.csv"

# Fixed threshold (HIGH STRINGENCY)
BRAIN_RPM_THRESHOLD = 10.0

# --------------------------------------------------
# Load shared miRNAs (TarBase ∩ miRTarBase)
# --------------------------------------------------
shared_df = pd.read_csv(shared_mirna_file)
shared_df.columns = shared_df.columns.str.strip()

shared_mirna_col = shared_df.columns[0]
shared_mirnas = set(shared_df[shared_mirna_col].astype(str).str.strip())

# --------------------------------------------------
# Load miRmine expression data
# --------------------------------------------------
mirmine = pd.read_csv(mirmine_file)
mirmine.columns = mirmine.columns.str.strip()

# miRNA identifier column (from your file)
mirna_col = "Mature miRNA ID"

# All brain RNA-seq sample columns
brain_cols = [c for c in mirmine.columns if "(Brain)" in c]

print(f"Detected {len(brain_cols)} brain samples")

# --------------------------------------------------
# Compute mean brain expression
# --------------------------------------------------
mirmine["Mean_Brain_RPM"] = mirmine[brain_cols].mean(axis=1)

# --------------------------------------------------
# Filter: shared miRNAs with expression ≥ 10 RPM
# --------------------------------------------------
brain_mirnas = mirmine[
    (mirmine[mirna_col].isin(shared_mirnas)) &
    (mirmine["Mean_Brain_RPM"] >= BRAIN_RPM_THRESHOLD)
]

# --------------------------------------------------
# Save outputs
# --------------------------------------------------
brain_mirnas[[mirna_col, "Mean_Brain_RPM"]].rename(
    columns={mirna_col: "miRNA"}
).to_csv(
    "Shared_BrainExpressed_miRNAs_RPM10.csv",
    index=False
)

pd.Series(
    sorted(brain_mirnas[mirna_col].unique()),
    name="miRNA"
).to_csv(
    "Shared_BrainExpressed_miRNA_list_RPM10.txt",
    index=False
)

# --------------------------------------------------
# Summary
# --------------------------------------------------
print("Brain-expression filtering completed.")
print(f"Shared miRNAs (input)                     : {len(shared_mirnas)}")
print(f"Brain-expressed miRNAs (mean RPM ≥ 10)    : {brain_mirnas[mirna_col].nunique()}")

