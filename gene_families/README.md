# Gene family analysis

Candidate gene family survey across *Mespilia globulus* ("Tuxedo") and five
comparator species (Hleu, Pliv, Spurp, Lvar, Lpic — same taxa as
`../synteny/`), built on top of an OrthoFinder orthology run: which families
were gained/lost/duplicated (CAFE), which show significant copy-number
shifts per branch, and how the candidates of interest are expressed across
*Mespilia* tissues.

Reorganised and annotated from `../Gene_Family_Analysis/` (left untouched
for provenance). Two exact-duplicate files present in the original
(`gene_family_raw_counts.csv` copied into both `Copy_number_heatmap/` and
`RNA_heatmap/`; `Base_branch_probabilities_tot.txt` copied into both `CAFE/`
and `Copy_number_heatmap/`) were deduplicated to one canonical copy each. A
typo'd directory name (`outut_small_families`) was corrected to
`output_small_families`.

**Caveat, as with `synteny/`:** paths were rewritten to the layout below and
not re-executed — smoke-test before relying on these scripts.

## Candidate family funnel

16 families were BLAST-searched against the 5 comparator species
(`data/blast_results/`): FMO, MSP130, NACHT, PKS, SM30, SM50, SOX, SRCR,
TLR, Wnt, fox, gcml, hox11:13, opsin, pax, rhol. Of those, 8 were carried
forward to curated multi-species FASTAs (`data/gene_family_fastas/`):
MSP130, PKS, SM30, SM50, SOX, Wnt, fox, opsin. Of those, 5 currently have
built gene trees (`data/gene_family_trees/`, plus a manually re-labelled
`annotated/` version used for figures): MSP130, SM30, SM50, SOX, fox.
SM30/SM50/MSP130 are biomineralization genes, SRCR/TLR/NACHT immune
receptors, SOX/fox/pax developmental transcription factors, PKS likely
relevant to blue/red pigment biosynthesis (see `../sex_determination/`).

**Note on sequence IDs:** FASTA and tree files use OrthoFinder's internal
`SpeciesN|...` sequence-ID convention (`Species1` … `Species5` for the
non-*Mespilia* taxa) rather than the species abbreviations used elsewhere in
this repo. Decoding `SpeciesN` → actual species requires
`OrthoFinder_output/WorkingDirectory/SpeciesIDs.txt`, which is not part of
this archive.

## Pipeline

1. **OrthoFinder** (external, not run from this repo) — orthology inference
   across all 6 species. Its raw output
   (`Orthogroups.tsv`, `Orthogroups.GeneCount.tsv`,
   `WorkingDirectory/SpeciesIDs.txt`) is **not included here** — two scripts
   below depend on it and can't currently be re-run (see their headers).
2. **CAFE** (external tool) — gene-family size evolution against the
   time-calibrated species tree (`data/species_tree/`), run twice: once on
   `data/cafe/small_families_filtered_cafe_input.tsv` and once on
   `large_families_filtered_cafe_input.tsv`, producing
   `data/cafe/output_small_families/` and `output_large_families/`
   respectively.
3. **`scripts/extract_sig_exporcontr.py`** — run once inside each CAFE
   output directory; filters families significant at p<0.01 and tallies
   per-branch expansions/contractions into
   `expansion_contraction_nochange_summary.csv`. The further-processed
   summary files sitting alongside it (`significant_expansion_contraction_summary.csv`,
   `total_sig_exp_contr_summary.csv`) look like a manual small+large merge
   step that isn't captured by any script here.
4. **`scripts/piechart_cafe.R`** — final figure: species tree with a pie
   chart per node/tip showing the significant-expansion/contraction/no-change
   split, from `total_sig_exp_contr_summary.csv`. Output: `figures/tree_node_pies.{pdf,svg}`.
5. **`scripts/plot_heatmap_with_sig.py`** — final figure: gene-family
   copy-number z-score heatmap (`data/copy_number/`) with CAFE-significance
   asterisks overlaid from `data/cafe/Base_branch_probabilities_tot.txt`.
   Output: `figures/gene_family_heatmap.svg`.
   **Possible bug, flagged not fixed:** the script's species→CAFE-column
   mapping has Lpic and Lvar swapped relative to the actual column order in
   `Base_branch_probabilities_tot.txt` — if real, the Lpic/Lvar asterisks in
   the rendered heatmap mark the wrong species. Check against the original
   CAFE run before trusting that figure.
6. **`scripts/tissue_heatmap.py`** — final figure: tissue-expression heatmap
   for the candidate families, from `data/expression/Mglob_gathered_tpms_tissues.txt`.
   Depends on the missing `Orthogroups.tsv` (see caveat above).
7. **`scripts/gene_copy_table.py`** — standalone, not part of the chain
   above: per-species singleton/duplicate/multicopy fractions. Also depends
   on a missing OrthoFinder file (`Orthogroups.GeneCount.tsv`).

`figures/` is currently empty — none of the above have been (re-)run since
reorganising.

## Directory layout

```
scripts/                        7 scripts, see pipeline above for run order/dependencies
data/
  species_tree/                 time-calibrated Newick tree (CAFE input)
  cafe/                         CAFE inputs + both small/large-family output dirs
  copy_number/                  per-species gene-family copy-number counts + z-scores
  expression/                   Mglob tissue TPM table
  blast_results/                {family}_vs_Species{1-5}.tsv, 16 families x 5 species
  gene_family_fastas/           curated multi-species FASTAs, 8 families
  gene_family_trees/            gene trees + annotated/ (relabelled for figures), 5 families
figures/                        empty until scripts are re-run
```
