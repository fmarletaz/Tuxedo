import pandas as pd

# === INPUT ===
input_file = "Orthogroups.GeneCount.tsv"  # adjust path if needed
output_file = "gene_copy_classification_fractions_present.tsv"

# === STEP 1: LOAD DATA ===
df = pd.read_csv(input_file, sep="\t", index_col=0)

# === STEP 2: CLASSIFY EACH ENTRY ===
def classify(count):
    if count == 0:
        return None  # species not present
    elif count == 1:
        return 'singleton'
    elif count == 2:
        return 'duplicate'
    else:
        return 'multicopy'

classified = df.applymap(classify)

# === STEP 3: CALCULATE FRACTIONS (ONLY WHERE SPECIES HAS GENE) ===
summary = {}
for species in df.columns:
    present = classified[species].dropna()
    total = len(present)
    counts = present.value_counts()
    summary[species] = {
        'singleton': counts.get('singleton', 0) / total,
        'duplicate': counts.get('duplicate', 0) / total,
        'multicopy': counts.get('multicopy', 0) / total
    }

# === STEP 4: SAVE TO FILE ===
summary_df = pd.DataFrame(summary).T
summary_df.index.name = 'Species'
summary_df.to_csv(output_file, sep='\t')

print(f"Filtered fraction table saved to {output_file}")
