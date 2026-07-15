#!/bin/bash

# Step 1: Extract clg mapping from Hleu file (user will provide this)
echo "Reading clg mapping from Hleu file..."
# Assuming the file will be called hleu_clg_mapping.txt with format: gene/chromosome \t clg
# User needs to provide this file

# Step 2: Create ortholog mappings (assuming ortholog files have gene1 \t gene2 format)
echo "Creating ortholog mappings..."

# Hleu to Mglob
awk '{print $1 "\t" $2}' mbh/Hleu_gene_Mglob_reci.txt > hleu_to_mglob.txt
awk '{print $2 "\t" $1}' mbh/Hleu_gene_Mglob_reci.txt > mglob_to_hleu.txt

# Mglob to Pliv
awk '{print $1 "\t" $2}' mbh/Mglob_Pliv_reci.txt > mglob_to_pliv.txt
awk '{print $2 "\t" $1}' mbh/Mglob_Pliv_reci.txt > pliv_to_mglob.txt

# Pliv to Spurp
awk '{print $1 "\t" $2}' mbh/Pliv_Spurp_reci.txt > pliv_to_spurp.txt
awk '{print $2 "\t" $1}' mbh/Pliv_Spurp_reci.txt > spurp_to_pliv.txt

# Spurp to Lpic
awk '{print $1 "\t" $2}' mbh/Spurp_Lpic_reci.txt > spurp_to_lpic.txt
awk '{print $2 "\t" $1}' mbh/Spurp_Lpic_reci.txt > lpic_to_spurp.txt

# Lpic to Lvar  
awk '{print $1 "\t" $2}' mbh/Lpic_Lvar_proteome_reci.txt > lpic_to_lvar.txt
awk '{print $2 "\t" $1}' mbh/Lpic_Lvar_proteome_reci.txt > lvar_to_lpic.txt

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
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' hleu_clg_mapping.txt beds/Holothuria_holleu_clean.bed > beds/Hleu_synt_with_clg.bed

# Add clg to mglob_output_syn.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' mglob_clg_mapping.txt beds/mglob_output_syn.bed > beds/mglob_output_syn_with_clg.bed

# Add clg to pliv_output.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' pliv_clg_mapping.txt beds/pliv_output.bed > beds/pliv_output_with_clg.bed

# Add clg to spur_synt_cl.bed (overwrite if it exists)
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' spurp_clg_mapping.txt beds/spur_synt_cl.bed > beds/spur_synt_cl_with_clg.bed

# Add clg to lpic_cleaned.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lpic_clg_mapping.txt beds/lpic_cleaned.bed > beds/lpic_cleaned_with_clg.bed

# Add clg to lvar_output.bed  
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lvar_clg_mapping.txt beds/lvar_output.bed > beds/lvar_output_with_clg.bed

# Step 5: Clean up temporary files
echo "Cleaning up temporary files..."
rm *_to_*.txt *_clg_mapping.txt

echo "Done! Files with clg column in beds/ directory:"
echo "  beds/Hleu_synt_with_clg.bed (source)"
echo "  beds/mglob_output_syn_with_clg.bed"
echo "  beds/pliv_output_with_clg.bed"
echo "  beds/spur_synt_cl_with_clg.bed"
echo "  beds/lpic_cleaned_with_clg.bed"
echo "  beds/lvar_output_with_clg.bed"
echo ""
echo "NOTE: You need to provide 'hleu_clg_mapping.txt' with format:"
echo "gene_id<TAB>clg"
echo "Example:"
echo "KAJ8017360<TAB>A"
echo "KAJ8017361<TAB>B"
