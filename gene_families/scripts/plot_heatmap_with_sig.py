"""
plot_heatmap_with_sig.py -- final figure, gene-family copy-number heatmap.

Draws a heatmap of gene-family copy-number z-scores across species
(../data/copy_number/gene_family_zscores.csv), annotated with an asterisk
wherever that family's orthogroup(s) show a CAFE-significant (p<0.01)
copy-number change on that species' branch
(../data/cafe/Base_branch_probabilities_tot.txt, the combined small+large
CAFE run -- see extract_sig_exporcontr.py's header for how the two runs
were combined). Output: ../figures/gene_family_heatmap.svg.

CAVEAT -- probable bug carried over from the original script: in
species_to_cafe_col below, Lpic and Lvar are mapped to each other's CAFE
column ("Lpic" -> "Lvar<5>", "Lvar" -> "Lpic<6>"), whereas
../data/cafe/Base_branch_probabilities_tot.txt's header order is
Hleu<1>, Mglob<2>, Pliv<3>, Spurp_prefixed<4>, Lvar<5>, Lpic<6>. If this
mapping is really swapped, the asterisks on the Lpic/Lvar columns of the
rendered heatmap mark the wrong species. Worth checking against the
original CAFE run before trusting that figure.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams["svg.fonttype"] = "none"

# === FILES ===
zscore_file = Path("../data/copy_number/gene_family_zscores.csv")
raw_file = Path("../data/copy_number/gene_family_raw_counts.csv")
cafe_file = Path("../data/cafe/Base_branch_probabilities_tot.txt")

# === CAFE column mapping (updated tree order) ===
# NOTE: see CAVEAT above -- Lpic/Lvar look swapped here relative to the
# actual column order in Base_branch_probabilities_tot.txt.
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
plt.savefig("../figures/gene_family_heatmap.svg", format="svg", bbox_inches="tight")
plt.show()
