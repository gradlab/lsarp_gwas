#! /usr/bin/env python
import os
import sys
import pandas as pd
import gzip

# unitig_to_itol.py ACGTAGTCATGC genename color

kmers = '/bulk/LSARP/genomics/analyses/GWAS/data/unitigs/s_aureus_unitigs/unitigs.txt'
kmer = sys.argv[1]
gene = sys.argv[2]
color = sys.argv[3]
with open(kmers, 'r') as dat:
    for line in dat:
        if kmer == line.rstrip().split('|')[0].strip():
            break
samples = line.rstrip().split('|')[1].split()
samples = [x.split(':')[0] for x in samples]

itol_header = '''
DATASET_COLORSTRIP
SEPARATOR COMMA
DATASET_LABEL,{0}
COLOR,{1}
LEGEND_TITLE,{0}
LEGEND_SHAPES,1
LEGEND_COLORS,{1}
LEGEND_LABELS,1
STRIP_WIDTH,25
BORDER_WIDTH,5
BORDER_COLOR,#CCCCCC
DATA
'''.format(gene,color)

with open('itol_{0}.txt'.format(gene), 'w') as out:
    out.write(itol_header)
    for sample in samples:
        out.write('{0},{1},1\n'.format(sample,color))
