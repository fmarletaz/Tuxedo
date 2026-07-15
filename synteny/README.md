# Synteny analysis

Pairwise macrosynteny analysis across six echinoderms, anchored on the focal
species of the manuscript, *Mespilia globulus* ("Tuxedo"), including a
within-species comparison of its two colour morphs (blue/male vs red/female,
linked to the sex-determining region described in `../sex_determination/`).

This folder is a reorganised, annotated version of the original working
directory `../Urchins_Synteny/` (left untouched for provenance). The
`local/` subfolder there was the actual final working directory — all
scripts `setwd()`'d into it — so its contents are what became `data/` and
`figures/` here. Everything outside `local/` in the original directory was
an earlier, superseded iteration of the same analysis and has been moved to
[`archive_superseded/`](archive_superseded/).

**Caveat:** paths inside the scripts were rewritten to match this layout and
data files were renamed for clarity, but the pipeline has not been
re-executed end to end after reorganising — do a smoke test (at least
`scripts/enrichment_and_dotplots.R` through `scripts/ideogram_figures.R` on
one species pair) before relying on it for the submission.

## Species abbreviations

| Code | Species | Role |
|---|---|---|
| Mglob | *Mespilia globulus* | Focal species ("Tuxedo"); `blue`/`red` = the two colour morphs |
| Hleu | *Holothuria leucospilota* | Outgroup (sea cucumber) |
| Pliv | *Paracentrotus lividus* | Comparator |
| Spur / Spurp | *Strongylocentrotus purpuratus* | Comparator |
| Lpic | *Lytechinus pictus* | Comparator |
| Lvar | *Lytechinus variegatus* | Comparator |

CLG = conserved/ancestral linkage group, the unit of macrosynteny used
throughout (analogous to bilaterian ALGs).

Species-tree topology (Newick, `data/species_tree/`):
`(Mglob,(Parliv,(Lytpic,Strpur)))` — not read by any script here, kept as a
supporting reference for the same four taxa.

## Pipeline

1. **`scripts/macrosynteny.py`** (library) + **`scripts/prep_synteny.py`**
   (CLI) — for one species pair, reduce single-copy reciprocal-best-hit
   orthologues (`data/orthologs/`) to genes on each species' largest
   chromosomes (`data/gene_beds/*_with_clg.bed`), rank them by chromosome
   and position, and write a raw synteny table to
   `data/synteny_tables/raw/{sp1}_vs_{sp2}_synt.txt`.
2. **`scripts/propagate_clg.sh`** (canonical) / **`propagate_clg_alt_spur_rooted.sh`**
   (alternate) — propagate CLG labels along the ortholog chain across all
   six taxa, producing the `*_with_clg.bed` files consumed by step 1.
   Documented for provenance; not directly re-runnable here (see header
   comments — most intermediate ortholog files were not preserved).
3. **`scripts/enrichment_and_dotplots.R`** (sourcing `synteny_functions.R`)
   — Fisher's-exact-test enrichment of each chromosome x CLG cell
   (Bonferroni-corrected), writing `data/synteny_tables/annotated/{sp1}_vs_{sp2}_syntr.txt`
   (adds an `scol` significance column) plus exploratory dotplots.
4. **`scripts/ideogram_figures.R`** — final RIdeogram figures (one SVG +
   PNG per species pair, colour-coded by CLG) in `figures/ideograms/`.

Run scripts from inside `scripts/`; all paths are relative to that
location.

## Directory layout

```
scripts/                        6 pipeline scripts, in run order above
data/
  species_tree/                 Newick tree for Mglob/Pliv/Lpic/Spur
  orthologs/                    reciprocal-best-hit files (incomplete, see caveat above)
  gene_beds/                    per-species gene BEDs, with CLG column where propagated
  synteny_tables/
    raw/                        stage-1 output ({sp1}_vs_{sp2}_synt.txt)
    annotated/                  stage-3a output ({sp1}_vs_{sp2}_syntr.txt), incl. one
                                 manually-patched variant (pliv_vs_spur_syntr_manualfix.txt,
                                 not currently wired into ideogram_figures.R)
    karyotype.txt
figures/
  dotplots_exploratory/         early per-pair dotplots (position- and rank-based)
  per_chromosome_pair/blue_vs_red_Mglob/   one PNG per best-matching scaffold/chromosome
                                            pair, blue vs red Mglob morphs
  ideograms/svg/, ideograms/png/           final manuscript-candidate figures
  ideograms/pdf_pairwise/                  earlier PDF dotplot versions of the same panels
archive_superseded/              earlier iterations of this analysis, kept for provenance:
  early_loose_synteny_tables/    top-level *_synt.txt files + ideogram_f.r (an older,
                                  multi-manuscript master script mixing in unrelated taxa)
  newlpict_to_spurp/
  prep_synteny_out_cucumber_to_tuxedo/
```

## Known loose ends (flagged, not resolved)

- `propagate_clg.sh` (Hleu-rooted) was designated canonical, but the two
  propagation scripts were not verified to reproduce the exact
  `data/gene_beds/*_with_clg.bed` checked in here — treat both as
  documentation of the intended logic rather than a verified recipe.
- `local/Hleu_synt_with_clg-spur_synt.txt` in the original directory was
  empty (0 bytes, an aborted run) and was not carried over.
- `data/synteny_tables/annotated/pliv_vs_spur_syntr_manualfix.txt` is an
  alternate, apparently hand-corrected version of `pliv_vs_spur_syntr.txt`
  that `ideogram_figures.R` does not currently read — check whether the fix
  should be the one used for the submitted figure.
