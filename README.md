# Tuxedo sea urchin (*Mespilia globulus*) data and supplement repository

Code and data supporting the manuscript *The tuxedo sea urchin Mespilia globulus: A fast-developing and tractable model for genomic and developmental biology*
Analyses are organised per top-level folder:

| Folder | What it covers | Status |
|---|---|---|
| [`sex_determination/`](sex_determination/) | Chr4 sequencing-coverage and heterozygosity comparison between blue-male and red-female morphs, identifying the sex-linked region | data + script present, no README yet |
| [`synteny/`](synteny/) | Pairwise macrosynteny across six echinoderms (outgroup *Holothuria leucospilota* through *Mespilia*, *Paracentrotus*, *Strongylocentrotus*, two *Lytechinus* spp.), plus the blue-vs-red *Mespilia* comparison | reorganized and annotated — see [`synteny/README.md`](synteny/README.md) for the full pipeline |
| [`gene_families/`](gene_families/) | Candidate gene family survey: BLAST search, curated gene trees, CAFE gene-family size evolution, copy-number and tissue-expression heatmaps | reorganized and annotated — see [`gene_families/README.md`](gene_families/README.md) for the full pipeline |

## Sex determination

Coverage and heterozygosity of 10kb windows across chr4 (and genome-wide) in
blue-male vs red-female individuals, comparing normalized male/female
coverage ratio (`sex2.R`) to identify the sex-linked region, plus
non-parametric tests (Kruskal-Wallis, Dunn's test, Wilcoxon effect size) for
which chromosome shows the strongest sex-biased signal. Output figures:
`Ratio_Chrom_*.pdf`, `chr4_cov_het.pdf`, `het_bias_*.pdf`. This is the same
blue/red *Mespilia* comparison also analysed at the macrosynteny level in
`synteny/` (`mglob_blue_vs_red_*`).

## Gene families

BLAST survey of 16 candidate gene families across *Mespilia* and 5
comparator species, narrowed to 8 curated gene-family FASTAs and 5 built
gene trees, plus CAFE gene-family-size evolution (run separately on small
and large families) and copy-number/tissue-expression heatmaps. See
[`gene_families/README.md`](gene_families/README.md) for the pipeline,
species-ID caveats, and a flagged likely bug (Lpic/Lvar column swap) in one
of the heatmap scripts.

## Provenance folders

- [`Urchins_Synteny/`](Urchins_Synteny/) — the original, unreorganized
  synteny working directory that `synteny/` was built from. Left untouched
  for provenance; not needed for the manuscript itself.
- [`Gene_Family_Analysis/`](Gene_Family_Analysis/) — the original,
  unreorganized gene-family working directory that `gene_families/` was
  built from. Left untouched for provenance; not needed for the manuscript
  itself.
- `archive_superseded/` — earlier iterations of the synteny pipeline,
  currently misplaced at the repo root (should live under
  `synteny/archive_superseded/`, see that folder's note) pending the same
  sync cleanup mentioned earlier.
- `sex_determination/Gene_Family_Analysis/` — a stray, now-redundant
  fragment left over from the earlier Dropbox sync issue (its 2 files are a
  subset of what's now in `Gene_Family_Analysis/`/`gene_families/`). Safe to
  delete once you've confirmed nothing else depends on it — I didn't remove
  it myself since it's not mine to discard.
