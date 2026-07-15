#!/bin/bash

# Step 1: Extract clg mapping from spur bed file (assuming spur_synt_cl.bed has clg in column 5)
echo "Extracting clg from spur bed file..."
awk '{print $4 "\t" $5}' beds/spur_synt_cl.bed > spur_clg_mapping.txt

# Step 2: Create ortholog mappings (assuming ortholog files have gene1 \t gene2 format)
echo "Creating ortholog mappings..."

# Spurp to Lpic
awk '{print $1 "\t" $2}' mbh/Spurp_Lpic_reci.txt > spurp_to_lpic.txt
awk '{print $2 "\t" $1}' mbh/Spurp_Lpic_reci.txt > lpic_to_spurp.txt

# Lpic to Lvar  
awk '{print $1 "\t" $2}' mbh/Lpic_Lvar_proteome_reci.txt > lpic_to_lvar.txt
awk '{print $2 "\t" $1}' mbh/Lpic_Lvar_proteome_reci.txt > lvar_to_lpic.txt

# Pliv to Spurp
awk '{print $1 "\t" $2}' mbh/Pliv_Spurp_reci.txt > pliv_to_spurp.txt
awk '{print $2 "\t" $1}' mbh/Pliv_Spurp_reci.txt > spurp_to_pliv.txt

# Mglob to Pliv
awk '{print $1 "\t" $2}' mbh/Mglob_Pliv_reci.txt > mglob_to_pliv.txt
awk '{print $2 "\t" $1}' mbh/Mglob_Pliv_reci.txt > pliv_to_mglob.txt

# Hleu to Mglob
awk '{print $1 "\t" $2}' mbh/Hleu_protein_Mglob_reci.txt > hleu_to_mglob.txt
awk '{print $2 "\t" $1}' mbh/Hleu_protein_Mglob_reci.txt > mglob_to_hleu.txt

# Step 3: Propagate clg forward (Spurp -> Lpic -> Lvar)
echo "Propagating clg forward..."

# Spurp to Lpic
awk 'NR==FNR {clg[$1]=$2; next} 
     {if ($1 in clg) print $2 "\t" clg[$1]}' spur_clg_mapping.txt spurp_to_lpic.txt > lpic_clg_mapping.txt

# Lpic to Lvar  
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' lpic_clg_mapping.txt lpic_to_lvar.txt > lvar_clg_mapping.txt

# Step 4: Propagate clg backward (Spurp -> Pliv -> Mglob -> Hleu)
echo "Propagating clg backward..."

# Spurp to Pliv
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' spur_clg_mapping.txt spurp_to_pliv.txt > pliv_clg_mapping.txt

# Pliv to Mglob
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' pliv_clg_mapping.txt pliv_to_mglob.txt > mglob_clg_mapping.txt

# Mglob to Hleu
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($1 in clg) print $2 "\t" clg[$1]}' mglob_clg_mapping.txt mglob_to_hleu.txt > hleu_clg_mapping.txt

# Step 5: Add clg column to each bed file
echo "Adding clg columns to bed files..."

# Add clg to lpic_cleaned.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lpic_clg_mapping.txt beds/lpic_cleaned.bed > beds/lpic_cleaned_with_clg.bed

# Add clg to lvar_output.bed  
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' lvar_clg_mapping.txt beds/lvar_output.bed > beds/lvar_output_with_clg.bed

# Add clg to pliv_output.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' pliv_clg_mapping.txt beds/pliv_output.bed > beds/pliv_output_with_clg.bed

# Add clg to mglob_output_syn.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' mglob_clg_mapping.txt beds/mglob_output_syn.bed > beds/mglob_output_syn_with_clg.bed

# Add clg to Hleu_synt.bed
awk 'NR==FNR {clg[$1]=$2; next}
     {if ($4 in clg) print $0 "\t" clg[$4]; else print $0 "\tNA"}' hleu_clg_mapping.txt beds/Hleu_synt.bed > beds/Hleu_synt_with_clg.bed

# Step 6: Clean up temporary files
echo "Cleaning up temporary files..."
rm *_to_*.txt *_clg_mapping.txt

echo "Done! Files with clg column in beds/ directory:"
echo "  beds/lpic_cleaned_with_clg.bed"
echo "  beds/lvar_output_with_clg.bed" 
echo "  beds/pliv_output_with_clg.bed"
echo "  beds/mglob_output_syn_with_clg.bed"
echo "  beds/Hleu_synt_with_clg.bed"
echo "  beds/spur_synt_cl.bed (original with clg)"
