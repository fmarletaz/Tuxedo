"""
extract_sig_exporcontr.py -- pipeline stage 2 (of 3), CAFE post-processing.

Run from inside a CAFE output directory (../data/cafe/output_small_families/
or ../data/cafe/output_large_families/, one run per family-size bin -- CAFE
was run separately on small and large gene families). Loads that run's
Base_family_results.txt (per-family p-value + significance call) and
Base_change.tab (per-branch copy-number change for every family), keeps only
families significant at p<0.01, and tallies per-branch counts of
expansions/contractions/no-change into expansion_contraction_nochange_summary.csv.

CAVEAT: the extra summary files sitting alongside this script's output
(significant_expansion_contraction_summary.csv, total_sig_exp_contr_summary.csv,
and the large_-prefixed equivalents) were not reproduced by this script as
checked in here -- they look like a further manual combination/filtering
step (e.g. merging the small+large runs) that wasn't captured. total_sig_exp_contr_summary.csv
is what piechart_cafe.R actually reads.
"""

import pandas as pd

# === STEP 1: Load and filter significant gene families ===
print("Loading family results...")
fam = pd.read_csv("Base_family_results.txt", sep="\t")

# Filter for p < 0.01 and significance == 'y'
sig_fams = fam[(fam.iloc[:, 1] < 0.01) & (fam.iloc[:, 2] == 'y')]
sig_ids = sig_fams.iloc[:, 0].tolist()
print(f"Found {len(sig_ids)} significant gene families.")

# === STEP 2: Load Base_change.tab and extract numeric values ===
print("Loading gene copy number changes...")
change = pd.read_csv("Base_change.tab", sep="\t")
change.rename(columns={change.columns[0]: "FamilyID"}, inplace=True)

# Numeric matrix for all families
numeric = change.drop(columns="FamilyID").apply(pd.to_numeric, errors="coerce")

# Subset only significant families
sig_change = change[change["FamilyID"].isin(sig_ids)]
sig_numeric = sig_change.drop(columns="FamilyID").apply(pd.to_numeric, errors="coerce")

# === STEP 3: Count expansions and contractions ===
print("Counting expansions and contractions...")
expansions = (sig_numeric > 0).sum()
contractions = (sig_numeric < 0).sum()

# === STEP 4: Count all families present at each node (non-zero change) ===
print("Counting families present at each node...")
total_present = (numeric != 0).sum()

# === STEP 5: Compute 'No Change' as remainder ===
print("Calculating no-change counts...")
no_change = total_present - expansions - contractions

# === STEP 6: Combine and save ===
summary = pd.DataFrame({
    "Node": expansions.index,
    "Significant_Expansions": expansions.values,
    "Significant_Contractions": contractions.values,
    "No_Change": no_change.values
})

output_file = "expansion_contraction_nochange_summary.csv"
summary.to_csv(output_file, index=False)
print(f"\nDone! Output saved to {output_file}")
