# ideogram_figures.R -- pipeline stage 3b (of 3), final figure generation.
#
# For each species pair, reads the enrichment-annotated synteny table
# (data/synteny_tables/annotated/*_syntr.txt, produced by
# enrichment_and_dotplots.R), keeps only significantly-enriched links
# (scol=='S'), builds an RIdeogram karyotype + synteny table coloured by
# CLG (conserved/ancestral linkage group), and writes the ideogram figure
# to figures/ideograms/svg/ (+ a PNG conversion in figures/ideograms/png/).
#
# Chromosome display order per species (k_*_ord vectors) was set by hand
# to align homologous chromosomes/CLGs visually across panels -- these
# orderings are manuscript-specific and were derived by inspecting the
# dotplots from enrichment_and_dotplots.R.
#
# Run from within scripts/ (paths below are relative to this file's location).

require(RIdeogram)
require(tidyverse)
library("wesanderson")


`%notin%` <- Negate(`%in%`)

ANNOT_DIR <- "../data/synteny_tables/annotated"
SVG_DIR <- "../figures/ideograms/svg"
PNG_DIR <- "../figures/ideograms/png"
dir.create(SVG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PNG_DIR, recursive = TRUE, showWarnings = FALSE)

# CLG (conserved/ancestral linkage group) colour palette, shared with
# enrichment_and_dotplots.R -- keep the two in sync if edited.
colors=data.frame(row.names = c('A1','A2','B1','B2+C2','B3','C1','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P', 'Q', 'R'),
                  col=c("A6CEE3","1F78B4","1B9E77","33A02C","B2DF8A","FB9A99","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","CAB2D6","E6AB02","A6761D","8DD3C7","666666","FF7F00","FDBF6F","BEBADA","ffd92f"))

# helper: after RIdeogram's convertSVG() writes a PNG next to the SVG, move
# it into PNG_DIR to keep the two formats separated as in the rest of the repo.
.moveConvertedPNG <- function(svg_path) {
  png_path <- sub("\\.svg$", ".png", svg_path)
  if (file.exists(png_path)) {
    file.rename(png_path, file.path(PNG_DIR, basename(png_path)))
  }
}

### ---- 1. Holothuria leucospilota (outgroup) vs M. globulus ----
lg=read.table(file.path(ANNOT_DIR, "hleu_vs_mglob_syntr.txt"),h=T)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='leucospilota',fill="leucospilota",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_hleu

lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='globulus',fill="globulus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_glob

k_hleu_ord=c("Hchr2","Hchr9","Hchr3", "Hchr22", "Hchr23","Hchr14",'Hchr17',
              "Hchr4", "Hchr5", "Hchr8","Hchr16", "Hchr15",
              "Hchr18", "Hchr21", "Hchr7",
             "Hchr6","Hchr20", "Hchr11","Hchr13", 'Hchr10',"Hchr12", 'Hchr1','Hchr19')

k_hleu=k_hleu[order(match(k_hleu$Chr,k_hleu_ord)),]

rbind(k_hleu,k_glob) -> chroms

lg$Species_1=match(lg$s1chr,k_hleu$Chr)
lg$Species_2=match(lg$s2chr,k_glob$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)

svg_out <- file.path(SVG_DIR, "hleu_vs_mglob_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 2. M. globulus vs Paracentrotus lividus ----
lg=read.table(file.path(ANNOT_DIR, "mglob_vs_pliv_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='globulus',fill="globulus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_mglob
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='lividus',fill="lividus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pliv

# Add chromosome ordering if needed
# k_mglob_ord=c("scaffold_1","scaffold_2","scaffold_3", ...) # Define based on your data
k_pliv_ord = c( "Scaffold_2100","Scaffold_174","Scaffold_2974","Scaffold_3425","Scaffold_3432","Scaffold_3434", "Scaffold_3435","Scaffold_3430","Scaffold_3432","Scaffold_218","Scaffold_3428",
                "Scaffold_2974","Scaffold_641", "Scaffold_3431","Scaffold_649", "Scaffold_1721", "Scaffold_3433",
                "Scaffold_3429",
                "Scaffold_674",
                "Scaffold_3426")
# k_mglob=k_mglob[order(match(k_mglob$Chr,k_mglob_ord)),]
k_pliv=k_pliv[order(match(k_pliv$Chr,k_pliv_ord)),]

rbind(k_mglob,k_pliv) -> chroms
lg$Species_1=match(lg$s1chr,k_mglob$Chr)
lg$Species_2=match(lg$s2chr,k_pliv$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)

svg_out <- file.path(SVG_DIR, "mglob_vs_pliv_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 3. P. lividus vs S. purpuratus ----
lg=read.table(file.path(ANNOT_DIR, "pliv_vs_spur_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='plividus',fill="plividus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pliv
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='purpuratus',fill="purpuratus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur

# Add chromosome ordering if needed
k_pliv_ord = c( "Scaffold_2100","Scaffold_174","Scaffold_2974","Scaffold_3425","Scaffold_3432","Scaffold_3434", "Scaffold_3435","Scaffold_3430","Scaffold_3432","Scaffold_218","Scaffold_3428",
                "Scaffold_2974","Scaffold_641", "Scaffold_3431","Scaffold_649", "Scaffold_1721", "Scaffold_3433",
                "Scaffold_3429",
                "Scaffold_674",
                "Scaffold_3426")
k_spur_ord=c("spur5_scaffold_6","spur5_scaffold_1",'spur5_scaffold_5',"spur5_scaffold_7",'spur5_scaffold_9','spur5_scaffold_15','spur5_scaffold_2',"spur5_scaffold_8","spur5_scaffold_14",'spur5_scaffold_21',
             'spur5_scaffold_11','spur5_scaffold_19','spur5_scaffold_16',
             'spur5_scaffold_4','spur5_scaffold_3','spur5_scaffold_10',"spur5_scaffold_20","spur5_scaffold_13","spur5_scaffold_12","spur5_scaffold_18",'spur5_scaffold_17') # Define based on your data
k_pliv=k_pliv[order(match(k_pliv$Chr,k_pliv_ord)),]
k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]

rbind(k_pliv,k_spur) -> chroms
lg$Species_1=match(lg$s1chr,k_pliv$Chr)
lg$Species_2=match(lg$s2chr,k_spur$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)

svg_out <- file.path(SVG_DIR, "pliv_vs_spur_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 4. S. purpuratus vs L. pictus ----
lg=read.table(file.path(ANNOT_DIR, "spur_vs_lpic_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='purpuratus',fill="purpuratus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='pictus',fill="pictus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lpic

k_lpic_ord=c("chr3","chr8","chr7", "chr1","scaffold_19", "chr5", "chr12",
             "chr2", "chr10", "chr4", "chr13", "chr14","chr15",
              "chr11", "chr18","chr6","chr17",
             "chr16", "chr9")
k_spur_ord=c("spur5_scaffold_6","spur5_scaffold_1",'spur5_scaffold_5',"spur5_scaffold_7",'spur5_scaffold_9','spur5_scaffold_15','spur5_scaffold_2',"spur5_scaffold_8","spur5_scaffold_14",'spur5_scaffold_21',
             'spur5_scaffold_11','spur5_scaffold_19','spur5_scaffold_16',
             'spur5_scaffold_4','spur5_scaffold_3','spur5_scaffold_10',"spur5_scaffold_20","spur5_scaffold_13","spur5_scaffold_12","spur5_scaffold_18",'spur5_scaffold_17') # Define based on your data
# Add spur chromosome ordering if needed

k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]
k_lpic=k_lpic[order(match(k_lpic$Chr,k_lpic_ord)),]

rbind(k_spur,k_lpic) -> chroms
lg$Species_1=match(lg$s1chr,k_spur$Chr)
lg$Species_2=match(lg$s2chr,k_lpic$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)
svg_out <- file.path(SVG_DIR, "spur_vs_lpic_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 5. L. pictus vs L. variegatus ----
lg=read.table(file.path(ANNOT_DIR, "lpic_vs_lvar_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='pictus',fill="pictus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lpic
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='variegatus',fill="variegatus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar

k_lpic_ord=c("chr3","chr8","chr7", "chr1","scaffold_19", "chr5", "chr12",
             "chr2", "chr10", "chr4", "chr13", "chr14","chr15",
             "chr11", "chr18","chr6","chr17",
             "chr16", "chr9")
k_lvar_ord=c("chr3","chr8","chr7", "chr1", "chr5", "chr12",
             "chr2", "chr10", "chr4", "chr13", "chr14",
             "chr15", "chr11", "chr18", "chr6",
             "chr17","chr16", "chr9", "chr19")

k_lpic=k_lpic[order(match(k_lpic$Chr,k_lpic_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]

rbind(k_lpic,k_lvar) -> chroms
lg$Species_1=match(lg$s1chr,k_lpic$Chr)
lg$Species_2=match(lg$s2chr,k_lvar$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)
svg_out <- file.path(SVG_DIR, "lpic_vs_lvar_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 6. S. purpuratus vs L. variegatus ----
lg=read.table(file.path(ANNOT_DIR, "spur_vs_lvar_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='purpuratus',fill="purpuratus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='variegatus',fill="variegatus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar

k_spur_ord=c("spur5_scaffold_6","spur5_scaffold_1",'spur5_scaffold_5',"spur5_scaffold_7",'spur5_scaffold_9','spur5_scaffold_15','spur5_scaffold_2',"spur5_scaffold_8","spur5_scaffold_14",'spur5_scaffold_21',
             'spur5_scaffold_11','spur5_scaffold_19','spur5_scaffold_16',
             'spur5_scaffold_4','spur5_scaffold_3','spur5_scaffold_10',"spur5_scaffold_20","spur5_scaffold_13","spur5_scaffold_12","spur5_scaffold_18",'spur5_scaffold_17') # Define based on your data
k_lvar_ord=c("chr3","chr8","chr7", "chr1", "chr5", "chr12",
             "chr2", "chr10", "chr4", "chr13", "chr14",
             "chr15", "chr11", "chr18", "chr6",
             "chr17","chr16", "chr9", "chr19")

k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]

rbind(k_spur,k_lvar) -> chroms
lg$Species_1=match(lg$s1chr,k_spur$Chr)
lg$Species_2=match(lg$s2chr,k_lvar$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)
svg_out <- file.path(SVG_DIR, "spur_vs_lvar_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 7. L. variegatus vs L. pictus (reciprocal direction, separate MBH run) ----
lg=read.table(file.path(ANNOT_DIR, "lvar_vs_lpic_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='variegatus',fill="variegatus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='pictus',fill="pictus",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lpic


k_lvar_ord=c("chr3","chr8","chr7", "chr1", "chr5", "chr12",
             "chr2", "chr10", "chr4", "chr13", "chr14",
             "chr15", "chr11", "chr18", "chr6",
             "chr17","chr16", "chr9", "chr19")

k_lpic=k_lpic[order(match(k_lpic$Chr,k_lpic_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]

rbind(k_lvar,k_lpic) -> chroms
lg$Species_1=match(lg$s1chr,k_lvar$Chr)
lg$Species_2=match(lg$s2chr,k_lpic$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)
svg_out <- file.path(SVG_DIR, "lvar_vs_lpic_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)

### ---- 8. Mespilia globulus: blue morph vs red morph ----
lg=read.table(file.path(ANNOT_DIR, "mglob_blue_vs_red_syntr.txt"),h=T)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>%
  mutate(Start=1,species='mglob blue',fill="mglob blue",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_bmglob
lg %>%  group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>%
  mutate(Start=1,species='mglob red',fill="mglob red",size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_rmglob

k_mglob_ord=c("scaffold_3","scaffold_1","scaffold_2", "scaffold_5", "scaffold_6", "scaffold_4",
                            "scaffold_7", "scaffold_9", "scaffold_8", "scaffold_11", "scaffold_14",
                            "scaffold_13", "scaffold_10", "scaffold_21", "scaffold_20",
                            "scaffold_12", "scaffold_17","scaffold_18", "scaffold_19", "scaffold_15", "scaffold_16")
k_bmglob=k_bmglob[order(match(k_bmglob$Chr,k_mglob_ord)),]
rbind(k_bmglob,k_rmglob) -> chroms
lg$Species_1=match(lg$s1chr,k_bmglob$Chr)
lg$Species_2=match(lg$s2chr,k_rmglob$Chr)
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt
lgt <- na.omit(lgt)
svg_out <- file.path(SVG_DIR, "mglob_blue_vs_red_syntr_ideogram.svg")
ideogram(karyotype = chroms, synteny = lgt, output = svg_out)
convertSVG(svg_out, device = "png")
.moveConvertedPNG(svg_out)
