from pathlib import Path
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================================================
# FILES
# =========================================================

raw_counts_file = Path("gene_family_raw_counts.csv")
orthogroups_file = Path("Orthogroups.tsv")
tpm_file = Path("Mglob_gathered_tpms_tissues.txt")

output_svg = Path("mglob_gene_family_expression_heatmap.svg")
output_matrix = Path("mglob_gene_family_expression_heatmap_matrix.csv")

# =========================================================
# SETTINGS
# =========================================================

gene_order = [
    "PKS", "FMO", "MSP130", "SM30", "SM50",
    "SOX", "fox", "pax",
    "NACHT", "SRCR", "TLR",
    "rhol", "opsin"
]

use_log2 = True
row_zscore = True

# Keep text editable in SVG
plt.rcParams["svg.fonttype"] = "none"

# =========================================================
# HELPERS
# =========================================================

def split_members(x):
    if pd.isna(x):
        return []
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return []
    return [i.strip() for i in re.split(r",\s*", x) if i.strip()]

def strip_isoform(tid):
    """
    Convert transcript IDs like MGLB00001.3 -> MGLB00001
    """
    if pd.isna(tid):
        return np.nan
    tid = str(tid).strip()
    m = re.match(r"^(MGLB\d+)(?:\.\d+)?$", tid)
    if m:
        return m.group(1)
    return np.nan

def rowwise_zscore(df):
    means = df.mean(axis=1)
    stds = df.std(axis=1).replace(0, np.nan)
    z = df.sub(means, axis=0).div(stds, axis=0)
    return z.fillna(0)

# =========================================================
# LOAD FAMILY -> ORTHOGROUP MAP
# =========================================================

raw_counts = pd.read_csv(raw_counts_file)

required_raw = {"Orthogroup", "GeneFamily"}
missing_raw = required_raw - set(raw_counts.columns)
if missing_raw:
    raise ValueError(f"Missing required columns in {raw_counts_file}: {missing_raw}")

wanted_ogs = (
    raw_counts.loc[raw_counts["GeneFamily"].isin(gene_order), ["Orthogroup", "GeneFamily"]]
    .drop_duplicates()
)

if wanted_ogs.empty:
    raise ValueError("No orthogroups found for the selected gene families.")

# =========================================================
# LOAD ORTHOGROUP MEMBERSHIP
# =========================================================

orthos = pd.read_csv(orthogroups_file, sep="\t")

required_orthos = {"Orthogroup", "Mglob"}
missing_orthos = required_orthos - set(orthos.columns)
if missing_orthos:
    raise ValueError(f"Missing required columns in {orthogroups_file}: {missing_orthos}")

orthos_sub = orthos.merge(wanted_ogs, on="Orthogroup", how="inner")

if orthos_sub.empty:
    raise ValueError("No overlap between Orthogroups.tsv and gene_family_raw_counts.csv.")

orthos_sub["Mglob_gene"] = orthos_sub["Mglob"].apply(split_members)
mglob_map = orthos_sub[["Orthogroup", "GeneFamily", "Mglob_gene"]].explode("Mglob_gene")
mglob_map = mglob_map.dropna(subset=["Mglob_gene"])
mglob_map = mglob_map[mglob_map["Mglob_gene"] != ""].drop_duplicates()

if mglob_map.empty:
    raise ValueError("No Mglob genes found after exploding Orthogroups.tsv.")

# =========================================================
# LOAD TPM TABLE
# =========================================================

tpm = pd.read_csv(tpm_file, sep="\t")

if "TID" not in tpm.columns:
    raise ValueError(f"{tpm_file} must contain a 'TID' column.")

tissue_cols = [c for c in tpm.columns if c != "TID"]
if not tissue_cols:
    raise ValueError(f"No tissue columns found in {tpm_file}.")

tpm["Mglob_gene"] = tpm["TID"].apply(strip_isoform)

if tpm["Mglob_gene"].isna().all():
    raise ValueError("Could not derive gene IDs from TPM transcript IDs.")

# =========================================================
# MERGE EXPRESSION WITH FAMILY MAP
# =========================================================

expr = mglob_map.merge(tpm, on="Mglob_gene", how="inner")

if expr.empty:
    print("Example Mglob gene IDs from Orthogroups.tsv:")
    print(mglob_map["Mglob_gene"].dropna().astype(str).head(10).tolist())
    print("Example transcript IDs from TPM table:")
    print(tpm["TID"].dropna().astype(str).head(10).tolist())
    print("Example derived gene IDs from TPM table:")
    print(tpm["Mglob_gene"].dropna().astype(str).head(10).tolist())
    raise ValueError("Merge returned 0 rows.")

expr = expr[["GeneFamily", "Mglob_gene", "TID"] + tissue_cols].drop_duplicates()

# =========================================================
# COLLAPSE TO GENE FAMILY LEVEL
# =========================================================
# First sum isoforms to gene level, then sum genes to family level

gene_expr = (
    expr.groupby(["GeneFamily", "Mglob_gene"], as_index=False)[tissue_cols]
    .sum()
)

family_expr = (
    gene_expr.groupby("GeneFamily", as_index=False)[tissue_cols]
    .sum()
)

family_expr["GeneFamily"] = pd.Categorical(
    family_expr["GeneFamily"],
    categories=gene_order,
    ordered=True
)

family_expr = family_expr.sort_values("GeneFamily")
matrix = family_expr.set_index("GeneFamily")[tissue_cols].copy()

# =========================================================
# TRANSFORM
# =========================================================

if use_log2:
    matrix = np.log2(matrix + 1)

if row_zscore:
    matrix = rowwise_zscore(matrix)
    cmap = "vlag"
    center = 0
    cbar_label = "Row z-score of log2(TPM + 1)"
else:
    cmap = "viridis"
    center = None
    cbar_label = "log2(TPM + 1)" if use_log2 else "TPM"

matrix.to_csv(output_matrix)

# =========================================================
# PLOT
# =========================================================

sns.set(font_scale=0.9)

plt.figure(figsize=(10, len(matrix) * 0.5))

ax = sns.heatmap(
    matrix,
    cmap=cmap,
    center=center,
    linewidths=0.5,
    cbar_kws={"label": cbar_label},
    square=False
)

ax.set_title("Mglob gene family tissue expression heatmap")
ax.set_xlabel("Tissue")
ax.set_ylabel("Gene family")

plt.tight_layout()
plt.savefig(output_svg, format="svg", bbox_inches="tight")
plt.show()

print(f"\nDone. Heatmap saved to: {output_svg}")
print(f"Matrix saved to: {output_matrix}")