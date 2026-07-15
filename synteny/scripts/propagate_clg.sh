#!/bin/bash
#
# propagate_clg.sh -- pipeline stage 1b (of 3), canonical CLG-propagation
# recipe (Hleu-rooted).
#
# Propagates CLG (conserved/ancestral linkage group) labels along the
# ortholog chain Hleu -> Mglob -> Pliv -> Spurp -> Lpic -> Lvar, using a
# hand-curated Hleu-to-CLG mapping as the root, and writes a `beds/*_with_clg.bed`
# for every species. These *_with_clg.bed files are the ones actually
# checked into ../data/gene_beds/ and consumed by prep_synteny.py.
#
# CAVEAT: this is the documented recipe for how ../data/gene_beds/*_with_clg.bed
# were produced, kept for provenance -- it is not directly re-runnable from
# this archive. Only two of the referenced ../data/orthologs/*_reci.txt
# files (Hleu_gene_Mglob_reci.txt, Hleu_protein_Spurp_reci.txt) and no
# hleu_clg_mapping.txt (the hand-curated root CLG assignment for Hleu) were
# preserved. See propagate_clg_alt_spur_rooted.sh for an alternate,
# Spur-rooted version of the same logic.
#
# Step 1: Extract clg mapping from Hleu file (user will provide this)
echo "Reading clg mapping from Hleu file..."
# Assuming the file will be called hleu_clg_mapping.txt with format: gene/chromosome \t clg
# User needs to provide this file

# Step 2: Create ortholog mappings (assuming ortholog files have gene1 \t gene2 format)
echo "Creating ortholog mappings..."

# Hleu to Mglob
awk '{print $1 "\t" $2}' ../data/orthologs/Hleu_gene_Mglob_reci.txt > hleu_to_mglob.txt
awk '{print $2 "\t" $1}' ../data/orthologs/Hleu_gene_Mglob_reci.txt > mglob_to_hleu.txt

# Mglob to Pliv
awk '{print $1 "\t" $2}' ../data/orthologs/Mglob_Pliv_reci.txt > mglob_to_pliv.txt
awk '{print $2 "\t" $1}' ../data/orthologs/Mglob_Pliv_reci.txt > pliv_to_mglob.txt

# Pliv to Spurp
awk '{print $1 "\t" $2}' ../data/orthologs/Pliv_Spurp_reci.txt > pliv_to_spurp.txt
awk '{print $2 "\t" $1}' ../data/orthologs/Pliv_Spurp_reci.txt > spurp_to_pliv.txt

# Spurp to Lpic
awk '{print $1 "\t" $2}' ../data/orthologs/Spurp_Lpic_reci.txt > spurp_to_lpic.txt
awk '{print $2 "\t" $1}' ../data/orthologs/Spurp_Lpic_reci.txt > lpic_to_spurp.txt

# Lpic to Lvar  
awk '{print $1 "\t" $2}' ../data/orthologs/Lpic_Lvar_proteome_reci.txt > lpic_to_lvar.txt
awk '{print $2 "\t" $1}' ../data/orthologs/Lpic_Lvar_proteome_reci.txt > lvar_to_lpic.txt

# Step 3: Propagate clg forward from Hleu through the chain
echo "Propagating clg from Hleu through the ortholog chain..."

# Hleu to Mglob
awk 'NR==FNR {clg[$1]=$2; next} 
     {if ($1 in clg) print $2 "\t" clg[$1]}' hleu_clg_mapping.txt hleu_to_mglob.txt > mglob_clg_mapping.txt

# Mglob to Pliv
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' mglob_clg_mapping.txt mglob_to_pliv.txt > pliv_clg_mapping.txt

# Pliv to Spurp
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' pliv_clg_mapping.txt pliv_to_spurp.txt > spurp_clg_mapping.txt

# Spurp to Lpic
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' spurp_clg_mapping.txt spurp_to_lpic.txt > lpic_clg_mapping.txt

# Lpic to Lvar
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' lpic_clg_mapping.txt lpic_to_lvar.txt > lvar_clg_mapping.txt

# Step 4: Add clg column to each bed file
echo "Adding clg columns to bed files..."

# Hleu already has clg (from original mapping file)
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' hleu_clg_mapping.txt ../data/gene_beds/Holothuria_holleu_clean.bed > ../data/gene_beds/Hleu_synt_with_clg.bed

# Add clg to mglob_output_syn.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' mglob_clg_mapping.txt ../data/gene_beds/mglob_output_syn.bed > ../data/gene_beds/mglob_output_syn_with_clg.bed

# Add clg to pliv_output.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' pliv_clg_mapping.txt ../data/gene_beds/pliv_output.bed > ../data/gene_beds/pliv_output_with_clg.bed

# Add clg to spur_synt_cl.bed (overwrite if it exists)
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' spurp_clg_mapping.txt ../data/gene_beds/spur_synt_cl.bed > ../data/gene_beds/spur_synt_cl_with_clg.bed

# Add clg to lpic_cleaned.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lpic_clg_mapping.txt ../data/gene_beds/lpic_cleaned.bed > ../data/gene_beds/lpic_cleaned_with_clg.bed

# Add clg to lvar_output.bed  
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lvar_clg_mapping.txt ../data/gene_beds/lvar_output.bed > ../data/gene_beds/lvar_output_with_clg.bed

# Step 5: Clean up temporary files
echo "Cleaning up temporary files..."
rm *_to_*.txt *_clg_mapping.txt

echo "Done! Files with clg column in ../data/gene_beds/ directory:"
echo "  ../data/gene_beds/Hleu_synt_with_clg.bed (source)"
echo "  ../data/gene_beds/mglob_output_syn_with_clg.bed"
echo "  ../data/gene_beds/pliv_output_with_clg.bed"
echo "  ../data/gene_beds/spur_synt_cl_with_clg.bed"
echo "  ../data/gene_beds/lpic_cleaned_with_clg.bed"
echo "  ../data/gene_beds/lvar_output_with_clg.bed"
echo ""
echo "NOTE: You need to provide 'hleu_clg_mapping.txt' with format:"
echo "gene_id<TAB>clg"
echo "Example:"
echo "KAJ8017360<TAB>A"
echo "KAJ8017361<TAB>B"
