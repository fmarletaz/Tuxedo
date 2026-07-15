# Tuxedo (*Mespilia globulus*) — manuscript supplement

Code and data supporting the manuscript on sex determination and genome
evolution in the sea urchin *Mespilia globulus*, informally "Tuxedo" for its
two colour morphs (blue/male, red/female). Three analyses, one per top-level
folder:

| Folder | What it covers | Status |
|---|---|---|
| [`sex_determination/`](sex_determination/) | Chr4 sequencing-coverage and heterozygosity comparison between blue-male and red-female morphs, identifying the sex-linked region | data + script present, no README yet |
| [`synteny/`](synteny/) | Pairwise macrosynteny across six echinoderms (outgroup *Holothuria leucospilota* through *Mespilia*, *Paracentrotus*, *Strongylocentrotus*, two *Lytechinus* spp.), plus the blue-vs-red *Mespilia* comparison | reorganized and annotated — see [`synteny/README.md`](synteny/README.md) for the full pipeline |
| `gene_families/` | OrthoFinder-based gene family analysis (orthogroups, gene trees, CAFE gene-family size evolution, RNA expression heatmaps) | **in progress**, see caveat below |

## Sex determination

Coverage and heterozygosity of 10kb windows across chr4 (and genome-wide) in
blue-male vs red-female individuals, comparing normalized male/female
coverage ratio (`sex2.R`) to identify the sex-linked region, plus
non-parametric tests (Kruskal-Wallis, Dunn's test, Wilcoxon effect size) for
which chromosome shows the strongest sex-biased signal. Output figures:
`Ratio_Chrom_*.pdf`, `chr4_cov_het.pdf`, `het_bias_*.pdf`. This is the same
blue/red *Mespilia* comparison also analysed at the macrosynteny level in
`synteny/` (`mglob_blue_vs_red_*`).

## Gene families — in progress

The `Gene_Family_Analysis` source data is still being synced/relocated (it
currently sits at `sex_determination/Gene_Family_Analysis/` — a stray
location left over from an in-progress Dropbox sync, not its intended home)
and only partially downloaded: of the OrthoFinder output, CAFE gene-family
evolution results, BLAST results, and analysis scripts, only the species
tree and an RNA-seq TPM table have synced so far. This folder will get the
same reorganize-and-annotate treatment as `synteny/` once the source data is
fully in place.

## Provenance folders

- [`Urchins_Synteny/`](Urchins_Synteny/) — the original, unreorganized
  synteny working directory that `synteny/` was built from. Left untouched
  for provenance; not needed for the manuscript itself.
- `archive_superseded/` — earlier iterations of the synteny pipeline,
  currently misplaced at the repo root (should live under
  `synteny/archive_superseded/`, see that folder's note) pending the same
  sync cleanup mentioned above.
