import os
import pandas as pd
from DeepPurpose import DTI as models
from DeepPurpose.dataset import load_broad_repurposing_hub

# ================== Configuration ==================
SAVE_PATH = './saved_path'
MORGAN_MODEL_NAME = 'Morgan_CNN_BindingDB'   # pretrained from DeepPurpose
MPNN_MODEL_NAME = 'mpnn_dti_model'           # your own trained model
RESULT_DIR = './result'

os.makedirs(SAVE_PATH, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

# ================== Helper: check model folders ==================
def ensure_model_folder(name, friendly_name=None):
    if friendly_name is None:
        friendly_name = name
    if not os.path.exists(name):
        print(f"❗ Model folder '{name}' not found in current directory.")
        print(f"   This is required for {friendly_name}.")
        print("   Make sure the folder exists here, e.g.:")
        print(f"   cp -r /path/to/{name} ./")
        exit(1)

# ================== Ensure model folders exist ==================
ensure_model_folder(MORGAN_MODEL_NAME, "Morgan_CNN_BindingDB pretrained model")
ensure_model_folder(MPNN_MODEL_NAME, "your MPNN DTI model")

# ================== Load models ==================
print("🔍 Loading pretrained Morgan_CNN_BindingDB model...")
net_morgan = models.model_pretrained(MORGAN_MODEL_NAME)
print("✅ Morgan_CNN_BindingDB model loaded!")

print("\n🔍 Loading your MPNN DTI model...")
net_mpnn = models.model_pretrained(MPNN_MODEL_NAME)
print("✅ mpnn_dti_model loaded!")

# ================== Load Broad Repurposing Hub ==================
print("\n📦 Loading Broad Repurposing Hub dataset...")
X_repurpose, drug_name, drug_cid = load_broad_repurposing_hub(SAVE_PATH)
print(f"Loaded {len(X_repurpose)} compounds.")

# ================== HUB GENE TARGET (EDIT THIS FOR EACH HUB) ==================
# Replace this block for each hub gene you want to analyze
target_seq = [
""">PHIP
MSCERKGLSELRSELYFLIARFLEDGPCQQAAQVLIREVAEKELLPRRTDWTGKEHPRTY
QNLVKYYRHLAPDHLLQICHRLGPLLEQEIPQSVPGVQTLLGAGRQSLLRTNKSCKHVVW
KGSALAALHCGRPPESPVNYGSPPSIDCNLHEADTLFSRKLNGKYRLERLVPTAVYQHMK
MHKRILGHLSSVYCVTFDRTGRRIFTGSDDCLVKIWATDDGRLLATLRGHAAEISDMAVN
YENTMIAAGSCDKMIRVWCLRTCAPLAVLQGHSASITSLQFSPLCSGSKRYLSSTGADGT
ICFWLWDAGTLKINPRPAKFTERPRPGVQMICSSFSAGGMFLATGSTDHIIRVYFFGSGQ
PEKISELEFHTDKVDSIQFSNTSNRFVSGSRDGTARIWQFKRREWKSILLDMATRPAGQN
LQGIEDKITKMKVTMVAWDRHDNTVITAVNNMTLKVWNSYTGQLIHVLMGHEDEVFVLEP
HPFDPRVLFSAGHDGNVIVWDLARGVKIRSYFNMVAGRKPIRSGLMELENRKQLSD
"""
]

target_name = "PHIP"  # change this per hub

# ================== Run Drug Repurposing: Morgan model ==================
print(f"\n===== Running Drug Repurposing for {target_name} with Morgan_CNN_BindingDB =====")
repurpose_morgan = models.repurpose(
    X_repurpose,
    target_seq,
    net_morgan,
    drug_name,
    target_name
)

# ================== Run Drug Repurposing: MPNN model ==================
print(f"\n===== Running Drug Repurposing for {target_name} with mpnn_dti_model =====")
repurpose_mpnn = models.repurpose(
    X_repurpose,
    target_seq,
    net_mpnn,
    drug_name,
    target_name
)

# ================== Build DataFrame with raw scores ==================
df = pd.DataFrame({
    'Drug': drug_name,
    'Drug_CID': drug_cid,
    f'{target_name}_Morgan_raw': repurpose_morgan,
    f'{target_name}_MPNN_raw': repurpose_mpnn,
})

# ================== Per-target min–max normalization ==================
def minmax_normalize(series):
    s_min = series.min()
    s_max = series.max()
    if s_max == s_min:
        return pd.Series([0.5] * len(series), index=series.index)
    return (series - s_min) / (s_max - s_min)

df[f'{target_name}_Morgan_norm'] = minmax_normalize(df[f'{target_name}_Morgan_raw'])
df[f'{target_name}_MPNN_norm']  = minmax_normalize(df[f'{target_name}_MPNN_raw'])

# Combined normalized score (simple mean of the two)
df[f'{target_name}_Combined_norm'] = (
    df[f'{target_name}_Morgan_norm'] + df[f'{target_name}_MPNN_norm']
) / 2.0

# ================== Rank drugs for THIS hub gene ==================
df_sorted_combined = df.sort_values(
    by=f'{target_name}_Combined_norm',
    ascending=False
)

df_sorted_morgan = df.sort_values(
    by=f'{target_name}_Morgan_raw',
    ascending=False
)

df_sorted_mpnn = df.sort_values(
    by=f'{target_name}_MPNN_raw',
    ascending=False
)

print(f"\n===== Top 5 drugs for {target_name} by Combined normalized score =====")
print(df_sorted_combined[['Drug', 'Drug_CID',
                          f'{target_name}_Morgan_norm',
                          f'{target_name}_MPNN_norm',
                          f'{target_name}_Combined_norm']].head(5))

# ================== Save results ==================
prefix = f"{target_name}_repurposing"

combined_file = os.path.join(RESULT_DIR, f"{prefix}_combined_norm.csv")
morgan_file   = os.path.join(RESULT_DIR, f"{prefix}_Morgan_raw_norm.csv")
mpnn_file     = os.path.join(RESULT_DIR, f"{prefix}_MPNN_raw_norm.csv")

df_sorted_combined.to_csv(combined_file, index=False)
df_sorted_morgan.to_csv(morgan_file, index=False)
df_sorted_mpnn.to_csv(mpnn_file, index=False)

print(f"\n💾 Saved combined (norm) ranking for {target_name}: {combined_file}")
print(f"💾 Saved Morgan (raw+norm) for {target_name}:   {morgan_file}")
print(f"💾 Saved MPNN (raw+norm) for {target_name}:     {mpnn_file}")
