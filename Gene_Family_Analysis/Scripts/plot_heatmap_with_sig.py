import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams["svg.fonttype"] = "none"

# === FILES ===
zscore_file = Path("Genes_of_interest/gene_family_zscores.csv")
raw_file = Path("Genes_of_interest/gene_family_raw_counts.csv")
cafe_file = Path("Cafe/Base_branch_probabilities_tot.txt")

# === CAFE column mapping (updated tree order) ===
species_to_cafe_col = {
    "Hleu": "Hleu<1>",
    "Mglob": "Mglob<2>",
    "Pliv": "Pliv<3>",
    "Spurp_prefixed": "Spurp_prefixed<4>",
    "Lpic": "Lvar<5>",
    "Lvar": "Lpic<6>"
}

species_order = list(species_to_cafe_col.keys())

# === Desired gene family order ===
gene_order = [
    "PKS", "FMO", "MSP130", "SM30", "SM50",
    "SOX", "fox", "pax",
    "NACHT", "SRCR", "TLR",
    "rhol", "opsin"
]

# === Load z-scores and raw counts ===
z_df = pd.read_csv(zscore_file).set_index("GeneFamily")
raw_df = pd.read_csv(raw_file)
orthogroup_map = raw_df[["Orthogroup", "GeneFamily"]]

# === Load CAFE significance table ===
cafe_df = pd.read_csv(cafe_file, sep="\t")

# === Collect significant orthogroups per species ===
sig_ogs_by_species = {
    col: set(cafe_df[cafe_df[col] < 0.01]["FamilyID"])
    for col in species_to_cafe_col.values()
}

# === Map GeneFamily → Orthogroups ===
gene_to_ogs = orthogroup_map.groupby("GeneFamily")["Orthogroup"].apply(set).to_dict()

# === Build annotation matrix ===
annotations = pd.DataFrame("", index=z_df.index, columns=z_df.columns)

for gene_family, ogs in gene_to_ogs.items():
    for species, cafe_col in species_to_cafe_col.items():
        if gene_family in annotations.index and cafe_col in sig_ogs_by_species:
            if len(ogs & sig_ogs_by_species[cafe_col]) > 0:
                annotations.loc[gene_family, species] = "*"

# === Reorder rows and columns ===
z_df = z_df.loc[gene_order, species_order]
annotations = annotations.loc[gene_order, species_order]

# === Plot ===
plt.figure(figsize=(10, len(z_df) * 0.5))
sns.set(font_scale=0.9)

ax = sns.heatmap(
    z_df,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    annot=annotations,
    fmt='',
    cbar_kws={"label": "Z-score (gene family copy number)"},
    square=False
)

plt.title("Gene Family Z-scores with CAFE Significance")
plt.xlabel("Species")
plt.ylabel("Gene Family")
plt.tight_layout()
plt.savefig("gene_family_heatmap.svg", format="svg", bbox_inches="tight")
plt.show()
