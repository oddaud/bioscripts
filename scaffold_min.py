#!/usr/bin/env python3

'''
Unfinished script meant to parse out all scaffolds > 1000 bp after annotation

Does not work for some files so far, i.e. .genes, trnascan, trnascan.fasta
and .16s
and .cmsearch

'''

# aed - 2018

from subprocess import call
import argparse
from Bio import SeqIO
import os
import pandas as pd


def get_under(fa_file, length=1000):
    under = set()
    with open(fa_file) as fa:
        for record in SeqIO.parse(fa, 'fasta'):
            if len(record) < length:
                under.add(str(record.id))
    return under

def clean_all(basename, underset):
    pass
    # Cleans all files in a directory

def clean_faa(faa, underset, f='faa'):
    ofile = faafile.split('.')[0]+'_scaffold_min_1000.fa.genes.'+f

    of = open(ofile, 'w')
    with open(faafile) as faa:
        for record in SeqIO.parse(faa, 'fasta'):
            rid = str(record.id)
            scaf = '_'.join(rid.split('_')[:-1])
            if scaf not in underset:
                of.write('>'+str(record.description)+'\n')
                of.write(str(record.seq)+'\n')

    of.close()


def clean_bsix(bsix, underset):
    df = pd.read_csv(bsix, header=None, sep='\t')
    ofile = bsix.split('.')[0]+'_scaffold_min_1000.fa.genes.faa-vs-uniprot.b6+'

    df.columns = ['gene', 'a', 'b', 'c', 'd', 'e',
                  'f', 'g', 'h', 'i', 'j', 'k', 'l']

    for index, row in df.iterrows():
        scaf = '_'.join(row['gene'].split('_')[:-1])
        if scaf in underset:
            df.drop(index, inplace=True)

    df.to_csv(ofile, index=False, header=None, sep='\t')


def main():
    fa = '/data9/jsantini2/analysis/audra_prevotella/almeida_assemblies/spain_SRR5579978/test.fa'

    b = '/data9/jsantini2/analysis/audra_prevotella/almeida_assemblies/spain_SRR5579978/test.b6+'
    uset = get_under(fa)
    print(len(uset))

    clean_bsix(b, uset)



if __name__ == '__main__':
    main()
