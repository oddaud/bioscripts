#!/usr/bin/env python3

# aed - 2018

from subprocess import call
import argparse
from Bio import SeqIO
import os


def run_prodigal(fasta, cmds, code=11, loud=False):
    cmd = ' -g '+str(code)+' '+cmds
    fnaout = cmd.split('-d')[1].split(' ')[1]
    cmdrun = 'prodigal -i '+fasta+cmd
    if loud:
        print('Running prodigal command: ', cmdrun)
    run_cmd(cmdrun)
    return fnaout


def calc_cd(fna, fasta):
    cds = 0
    seq_len = 0
    for record in SeqIO.parse(fna, 'fasta'):
        cds += len(record.seq)

    for record in SeqIO.parse(fasta, 'fasta'):
        seq_len += len(record.seq)

    return cds/seq_len

def write_tsv():
    pass

def run_cmd(cmd, dry=False, shell=True):
    if shell:
        #print(cmd)
        if not dry: call(cmd,shell=True)
    else:
        #print(' '.join(cmd))
        if not dry: call(cmd)
        return


def argparser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=False)

    # General Arguments
    GenArgs = parser.add_argument_group('GENERAL ARGUMENTS')
    GenArgs.add_argument('-h', action="help",
                         help="show this help message and exit")

    GenArgs.add_argument("-f", "--fasta", required=True,
                         help='path to phage fasta file')

    GenArgs.add_argument("-c", "--prodigal_cmds",
                         help='custom prodigal commands')

    GenArgs.add_argument("-t", "--tsv",
                         help='file to APPEND coding densities to')

    GenArgs.add_argument("-cd", "--coding_density",
                         help='Minimum coding density with code 11 to trigger \
                         check for alternate coding',
                         default=0.9)

    args = parser.parse_args()
    return args


def main():
    a = argparser()
    f = a.fasta
    t = a.tsv
    mincd = a.coding_density

    # Build the commands
    if a.prodigal_cmds is not None:
        # Future addition: check for sensicle cmds
        cmds = a.prodigal_cmds

    else:
        o = f + '.genes'
        d = f + '.genes.faa'
        a = f + '.genes.fna'
        cmds = '-o %s -a %s -d %s -m -p single -q' % (o, a, d)

    # Run prodigal
    fna = run_prodigal(f, cmds)
    cd11 = calc_cd(fna, f)

    if cd11 < mincd:
        # Candidate for alt coding
        c = '-o temp.genes -d temp.fna -a temp.faa -m -p single -q'
        recoded_fna = run_prodigal(f, c, 15)

        cd15 = calc_cd(recoded_fna, f)
        if cd11 < cd15:
            # Code 15 is correct, rename temp files to overwrite code 11 files
            os.rename('temp.genes', o)
            os.rename('temp.fna', d)
            os.rename('temp.faa', a)
            code = 15
        else:
            code = 11
    else:
        cd15 = 'NA'
        code = 11

    if t is not None:
        of = open(t, 'a')
        of.write(f+'\t'+str(cd11)+'\t'+str(cd15)+'\t'+str(code)+'\n')
    else:
        print('Working on ', f)
        print(f+'\t'+str(cd11)+'\t'+str(cd15)+'\t'+str(code))


if __name__ == '__main__':
    main()
