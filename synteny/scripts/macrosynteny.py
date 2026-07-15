#!/usr/bin/env python
"""
macrosynteny.py -- core library for pairwise macrosynteny analysis.

Given gene coordinates (from a GTF or a 5-column BED with an optional CLG
column) for two species and a list of single-copy reciprocal-best-hit
orthologues between them, this module reduces the ortholog list to genes
that lie on the largest chromosomes/scaffolds of each species, orders them
by chromosome and position, and reports each orthologous pair together with
its rank-position index (used downstream for dotplots/ideograms) and its
CLG (conserved/ancestral linkage group) label.

Not run directly -- imported by prep_synteny.py, which is the CLI driver
for a single species pair. Pipeline stage: 1 (of 3), see ../README.md.
"""

import sys, csv
from pybedtools import BedTool
from collections import defaultdict


def geneLocSE(gtf):
    """Gene -> (chrom, start, end, strand) from a GTF, spanning all exons."""
    exons = BedTool(gtf)
    genePos = defaultdict(list)
    geneChr = {}
    geneStr = {}
    for lg in exons:
        gn = lg.attrs['gene_id']
        if lg[2] == 'exon':
            genePos[gn].append(int(lg.start))
            genePos[gn].append(int(lg.end))
            geneStr[gn] = lg.strand
            geneChr[gn] = str(lg.chrom)
    geneLoc = {}
    for gene in genePos:
        geneLoc[gene] = (geneChr[gene], min(genePos[gene]), max(genePos[gene]), geneStr[gn])
    return geneLoc


def geneLoc(gtf):
    """Gene -> (5'-most position, chrom) from a GTF, strand-aware."""
    exons = BedTool(gtf)
    genePos = defaultdict(list)
    geneChr = {}
    geneStr = {}
    for lg in exons:
        gn = lg.attrs['gene_id']
        if lg[2] == 'exon':
            if lg.strand == '+': genePos[gn].append(int(lg.start))
            if lg.strand == '-': genePos[gn].append(int(lg.end))
            geneStr[gn] = lg.strand
            geneChr[gn] = str(lg.chrom)
    geneLoc = {}
    for gene in genePos:
        if geneStr[gene] == '+':
            geneLoc[gene] = (min(genePos[gene]), geneChr[gene])
        elif geneStr[gene] == '-':
            geneLoc[gene] = (max(genePos[gene]), geneChr[gene])
    return geneLoc


def chrSize(genpos):
    """chrom -> max gene position + 5kb padding, from a gene->(pos,chrom) map."""
    byChr = defaultdict(list)
    for gen in genpos:
        byChr[genpos[gen][1]].append(genpos[gen][0])
    chrSize = {}
    for chrm in byChr:
        chrSize[chrm] = max(byChr[chrm]) + 5000
    return chrSize


def geneSort(genpos, pr, chrS):
    """Rank a set of genes `pr` by (chromosome size desc, position asc); gene -> (gene,pos,chrom,rank)."""
    pr = set(pr)
    red = []
    for g in pr:
        if g in genpos:
            red.append((g, genpos[g][0], genpos[g][1]))
        else:
            print(g, "coordinates not found!!")
    red = sorted(red, key=lambda x: (-chrS[x[2]], x[1]))
    redN = dict((a[0], (a[0], a[1], a[2], i)) for i, a in enumerate(red))
    return redN

def syntComp(s1Pos,s2Pos,pairs,s1ChrK,s2ChrK):
    rpairs=[(g1,g2) for g1,g2 in pairs if g1 in s1Pos and g2 in s2Pos]
    print(len(rpairs),'pairs of genes!')
    s1ChrS=chrSize(s1Pos)
    s2ChrS=chrSize(s2Pos)
    s1RedN=geneSort(s1Pos,[p[0] for p in rpairs],s1ChrS)
    s2RedN=geneSort(s2Pos,[p[1] for p in rpairs],s2ChrS)
    mergN=[]
    for g1,g2 in rpairs:
        if g1 in s1RedN and g2 in s2RedN:
            if s1RedN[g1][2] in s1ChrK and s2RedN[g2][2] in s2ChrK:
                mergN.append(list(s1RedN[g1])+list(s2RedN[g2]))
    mergNS=sorted(mergN,key=lambda x: (-s1ChrS[x[2]],-s2ChrS[x[6]],x[3],x[7]))
    return mergNS

def selectChromSize(chromSize,nchrom):
    """Names of the `nchrom` largest chromosomes/scaffolds."""
    return [str(c[0]) for c in sorted(chromSize.items(),key=lambda x:x[1],reverse=True)[0:nchrom]]

def selectChromGenes(chromGenes,cutoff):
    """Names of chromosomes/scaffolds carrying more than `cutoff` genes."""
    sel=[str(c[0]) for c in sorted(chromGenes.items(),key=lambda x:x[1],reverse=True) if c[1]>cutoff]
    print(len(sel))
    return sel

def chrGenes(genepos):
    """chrom -> gene count, from a gene->(pos,chrom) map."""
    chrGenes=defaultdict(int)
    for g in genepos:
        pos,chrm=genepos[g]
        chrGenes[chrm]+=1
    return chrGenes

def printSynt(spp_synt,filename):
    """Write the synteny table (with CLG column if present) as TSV."""
    with open(filename,'w') as out:
        if len(spp_synt[0])==9:
            out.write("s1g\ts1gp\ts1chr\ts1gi\ts2g\ts2gp\ts2chr\ts2gi\tclg\n")
        else:
            out.write("s1g\ts1gp\ts1chr\ts1gi\ts2g\ts2gp\ts2chr\ts2gi\n")
        for p in spp_synt:
            out.write('\t'.join(map(str,p))+'\n')

def printSyntRev(spp_synt,filename):
    """Write the synteny table with species 1/2 columns swapped."""
    with open(filename,'w') as out:
        out.write("s1g\ts1gp\ts1chr\ts1gi\ts2g\ts2gp\ts2chr\ts2gi\tclg\n")
        for p in spp_synt:
            pr=p[4:8]+p[0:4]+[p[8]]
            out.write('\t'.join(map(str,pr))+'\n')

def loadBed(fileName):
    """Parse a 4- or 5-column BED (chrom,start,end,gene[,clg]) -> (gene->(pos,chrom), gene->clg)."""
    Pos,CLG={},{}
    for i,rc in enumerate(csv.reader(open(fileName),delimiter='\t')):
        if i==0: continue
        Pos[rc[3]]=(int(rc[1]),rc[0])
        if len(rc)==4:
            CLG[rc[3]]='NA'
        else:
            if not rc[4]=='NA':
                CLG[rc[3]]=rc[4]
    return Pos,CLG

def loadCoord(fileName):
    """Like loadBed but for a (gene,chrom,pos[,...,clg]) coordinate table instead of a BED."""
    Pos,CLG={},{}
    for i,rc in enumerate(csv.reader(open(fileName),delimiter='\t')):
        if i==0: continue
        Pos[rc[0]]=(int(rc[2]),rc[1])
        if len(rc)==4:
            CLG[rc[0]]='NA'
        else:
            if not rc[4]=='NA':
                CLG[rc[0]]=rc[4]
    return Pos,CLG

def checkChrm(chrSize,chrGenes,cmx=50):
    """Print the top `cmx` chromosomes/scaffolds by size, with their gene counts."""
    for i,(ch,sz) in enumerate(sorted(chrSize.items(),key=lambda x:x[1],reverse=True)):
        print(i,ch,sz,chrGenes[ch])
        if i==cmx: break

def loadPairs(orthfile,rev=0):
    """Load a 2-column reciprocal-best-hit ortholog file as (gene1,gene2) tuples; rev=1 swaps the columns."""
    pairs=[]
    for rc in csv.reader(open(orthfile),delimiter='\t'):
        if rev==0:
            pairs.append((rc[0],rc[1]))
        elif rev==1:
            pairs.append((rc[1],rc[0]))
    return pairs

def addCLG(synt, refclg1, refclg2):
    """
    Add CLG annotations from both species' BED files.
    Prioritizes non-NA values and combines when both species have CLG data.
    """
    clgSynt = synt
    for op in synt:
        gene1_id = op[0]  # species 1 gene
        gene2_id = op[4]  # species 2 gene

        # Get CLG values, defaulting to 'NA' if not found
        clg1 = refclg1.get(gene1_id, 'NA')
        clg2 = refclg2.get(gene2_id, 'NA')

        # Determine final CLG value based on priority rules
        if clg1 != 'NA' and clg2 != 'NA':
            # Both have CLG - combine them
            final_clg = f"{clg1}|{clg2}"
        elif clg1 != 'NA':
            # Only species 1 has CLG
            final_clg = clg1
        elif clg2 != 'NA':
            # Only species 2 has CLG
            final_clg = clg2
        else:
            # Neither has CLG
            final_clg = 'NA'

        op.append(final_clg)

    return clgSynt
