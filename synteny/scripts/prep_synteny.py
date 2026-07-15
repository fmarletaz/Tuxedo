#!/usr/bin/env python
"""
prep_synteny.py -- CLI driver for a single pairwise macrosynteny comparison.

Reads two 4/5-column BED files (chrom, start, end, gene[, CLG]) and a
2-column reciprocal-best-hit ortholog file, restricts to genes on the N
largest chromosomes/scaffolds of each species, and writes a synteny table
"{bed1 prefix}-{bed2 prefix}_synt.txt" (gene, position, chromosome, rank
index for each species, plus CLG) into the current working directory.

Pipeline stage: 2 (of 3) -- consumes data/gene_beds/*.bed and
data/orthologs/*.txt, produces the "raw" tables in
data/synteny_tables/raw/. Requires the CLG column already present in the
input BEDs (see propagate_clg.sh, stage 1b). Example:

    python prep_synteny.py -b1 ../data/gene_beds/Hleu_synt_with_clg.bed \\
        -b2 ../data/gene_beds/mglob_output_syn_with_clg.bed \\
        -n1 23 -n2 21 -m ../data/orthologs/Hleu_gene_Mglob_reci.txt
"""

import sys
import argparse
from os import path

from macrosynteny import geneLoc
from macrosynteny import loadCoord
from macrosynteny import loadBed
from macrosynteny import loadPairs
from macrosynteny import chrSize
from macrosynteny import syntComp
from macrosynteny import addCLG
from macrosynteny import printSynt
from macrosynteny import selectChromSize


parser = argparse.ArgumentParser(description="Prepared a synteny file for a pair of species")
# Optional arguments with flags
parser.add_argument('-b1', '--bed1', help='Bed woth Gene location in species 1')
parser.add_argument('-b2', '--bed2', help='Bed woth Gene location in species 2')
parser.add_argument('-n1', '--nchr1', type=int, help='Number of chromosomes in species 1')
parser.add_argument('-n2', '--nchr2', type=int, help='Number of chromosomes in species 2')
parser.add_argument('-m', '--mbh', help='Single-copy orthologues between species 1 and 2')

#parser.add_argument('-o', '--output', help='Output file path')
#parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
#parser.add_argument('-n', '--number', type=int, default=1, help='Number of times to repeat')

args = parser.parse_args()


#species 1
n_chr1=args.nchr1
print(n_chr1)
pfx1=path.basename(args.bed1).split('.')[0]
Pos1,ALG1=loadBed(args.bed1)
chr1=selectChromSize(chrSize(Pos1),n_chr1)
print('Sp1',list(Pos1.items())[0:5])
print('Chr1',chr1)
# species 2
n_chr2=args.nchr2
pfx2=path.basename(args.bed2).split('.')[0]
Pos2,ALG2=loadBed(args.bed2)
chr2=selectChromSize(chrSize(Pos2),n_chr2)
print('Sp2',list(Pos2.items())[0:5],chr2[0:5],'...')
print('Chr2',chr2)
# load synteny comparisons
mbh=loadPairs(args.mbh) 
print('Pairs:',len(mbh),mbh[0:5])
# compute synteny 
synt=syntComp(Pos1,Pos2,mbh,chr1,chr2)
addCLG(synt, ALG1, ALG2)
#print(synt[0:5])
printSynt(synt,f"{pfx1}-{pfx2}_synt.txt")
