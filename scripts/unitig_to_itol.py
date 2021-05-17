#! /usr/bin/env python
import os
import sys
import pandas as pd
import gzip
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Create ITOL annotation file for a unitig")
    parser.add_argument("unitig_file", help="Unitig presence/absence file")
    parser.add_argument("unitig", help="Specific unitig of interest")
    parser.add_argument("unitig_name", help="Unitig name (e.g. what gene does it map to?)")
    parser.add_argument("color", help="HEX code for color of annotation")
    return parser.parse_args()

args = get_args()

kmers = args.unitig_file
kmer = args.unitig
gene = args.unitig_name
color = args.color

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
