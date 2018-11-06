'''
Script to annotate phage genes based on NR database, and then alignment
with homologs and HHPRED search

Created by Audra Devoto on March 24. 2018

Project: Mega Phage
'''
import numpy as np
import pandas as pd
from Bio import SeqIO
import pickle, os, sys, logging
from Bio.Alphabet import generic_dna
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict
import subprocess
from subprocess import call
import argparse

def generate_gbk(faa, ffasta, name):
    '''
    Given a fasta file and it's corresponding gene predictions, write a genbank file.
    name = a unique identifier < 16 characters long
    All genes are labeled as unknown. Returns that output file
    '''
    outfile = faa.split('.')[0]+'.gbk'
    faa_handle = open(faa, 'r')
    fasta_handle = open(ffasta, 'r')
    out_handle = open(outfile, 'w')

    for frecord in SeqIO.parse(fasta_handle, 'fasta'):
        sequence_string = str(frecord.seq)
        sequence_object = Seq(sequence_string, IUPAC.protein)

        newrecord = SeqRecord(sequence_object, id = str(frecord.id), name = name)

        for record in SeqIO.parse(faa_handle, 'fasta'):
            name = str(record.id)
            desc = str(record.description)
            start = int(desc.split('#')[1].rstrip())
            end = int(desc.split('#')[2].rstrip())
            strand = int(desc.split('#')[3])
            fname = OrderedDict([('locus_tag', name)])
            feature = SeqFeature(FeatureLocation(start=start, end=end),
                                                strand = strand, id = 'unknown',
                                                type='gene', qualifiers=fname)
            newrecord.features.append(feature)

        SeqIO.write(newrecord, out_handle, 'genbank')

    faa_handle.close()
    fasta_handle.close()
    out_handle.close()
    return outfile

def run_prodigal(fasta):
    faaout = fasta.split('.')[0]+'.faa'
    if os.path.isfile(faaout):
        print('\tPredicted genes found. Running analyses on \'%s\'' % (faaout,))
    else:
        ofile = fasta.split('.')[0]+'.gb'

        cmd = 'prodigal -i '+fasta+' -a '+faaout+' -m -o '+ofile

        print('Running command: ', cmd)
        run_cmd(cmd)
        os.remove(ofile)
        print('\n Prodigal run complete')
    return faaout

#WRITTEN BY MATT OLM
def run_cmd(cmd,dry=False,shell=True):
	if shell:
		print(cmd)
		if not dry: call(cmd,shell=True)
	else:
		print(' '.join(cmd))
		if not dry: call(cmd)
	return


def fix_psiblast(wd, filename):
    '''
    blast names often have commas in them. This function 'fixes' a psiblast
    output and replaces the commas in blast names with '--'
    '''
    of_name = wd.rstrip('/')+'/'+filename.split('/')[-1].rstrip('.csv')+'.fixed.csv'
    if os.path.isfile(of_name):
        print('Fixed version found! proceeding')
        return of_name

    of = open(of_name, 'w')
    with open(filename, 'r') as psi:
        for line in psi:
            if len(line.split(','))>8:
                last = '--'.join(line.split(',')[7:])

                newline = ','.join(line.split(',')[:7])+','+last
                of.write(newline)
            else:
                of.write(line)
    of.close()
    return of_name

def parse_psiblast(csv, evalue, pident):
    '''
    Takes as input a csv file (in the right format)
    returns a pandas dataframe with additional columns and filtered to
    specified evalue and percent identity of hits.
    '''
    evalue = float(evalue)
    pident = float(evalue)

    with open(csv, 'r') as f:
        df = pd.read_csv(csv, sep =',', header = None)

        df.columns = ['qseqid','sseqid','evalue','bitscore','length'
                    ,'pident','qcovs','subj_name']

        df = df[df.qseqid != 'Search has CONVERGED!']

        df_filtered = df[df['evalue']<evalue]
        df_filtered=df_filtered[df_filtered['pident']>pident]

        df_filtered['phage']=['_'.join(p.split('_')[:-1]) for p in df_filtered['qseqid']]
        df_filtered['hit']=[p.split('[')[0].rstrip() for p in df_filtered['subj_name']]

        return df_filtered

def choose_hit(df):
    '''
    Given a dataframe containing records for lots of phage:
    for each phage, chose the best hit for each feature. This is simply
    done based on sorting by evalue and pident (in the case of an evalue tie)

    Concat all the winners into one dataframe and return it.
    This data contains the winner for each feature from every phage.

    In other words, the ['qseqid'] column is nonredundant.
    '''

    phage = set(list(df['phage']))
    winners = []
    for p in phage:
        dfp = df[df['phage']==p]
        features = set(list(dfp['qseqid']))
        for f in features:
            dff = dfp[dfp['qseqid']==f]
            dff = dff.sort_values(['evalue', 'pident'], axis=0, ascending=True, inplace=False, kind='quicksort')
            winner = dff[:1]
            winners.append(winner)
    final_winners = pd.concat(winners)
    return final_winners

def sort_unknowns(df):
    #Silences a warning from pandas
    pd.options.mode.chained_assignment = None
    df.fillna(0)
    unk_terms = ['hypothetical','PREDICTED:' ,'uncharacterized','Uncharacterised',
    'unknown', 'unnamed', 'predicted protein', 'putative phage related protein',
    'Phage conserved protein']

    unknowns = df[df['subj_name'].str.contains('|'.join(unk_terms), case = False, na = False)]
    unknowns.loc[:,'funx'] = str('unknown')


    knowns = df[~df['subj_name'].str.contains('|'.join(unk_terms), case = False, na = False)]
    knowns.loc[:,'funx'] = str('known')

    return pd.concat([knowns, unknowns])

def return_feature_dict(knowns):
    if len(list(knowns['qseqid'])) == len(set(knowns['qseqid'])):
        knowns_dict = knowns.set_index('qseqid')['hit'].to_dict()
    else:
        print('uh oh, you must select the winners first')
    return knowns_dict

def add_NR_annotations(gb, features_dict, ofile=None):
    '''
    parse the genbank file
    If there is a new annotation that needs to be added, write it!
    Add new features on top of rather than replacing the original "unknowns"
    '''
    if ofile == None:
        ofile =gb.split('.')[0]+'.nr_added.gbk'

    count = 0
    for seq_record in SeqIO.parse(open(gb,"r"), "genbank"):
        count+=1
        for f in seq_record.features:
            if f.type == 'gene':
                fname = f.qualifiers['locus_tag'][0]
                s = f.location.start
                e = f.location.end
                strand = f.strand

                if fname in features_dict:
                    product = features_dict[fname]
                    quals = OrderedDict([('locus_tag', fname), ('product', product)])

                    newfeature = SeqFeature(FeatureLocation(start=s, end=e),
                                            strand = strand,
                                            type='NR',
                                            qualifiers=quals)

                    seq_record.features.append(newfeature)


        SeqIO.write(seq_record, ofile, 'genbank')
    if count >1:
        print('Oh no, it\'s a multi genbank...')

    return ofile

def peek(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return first

def check_input(infile):
    if os.path.isfile(infile):
        if infile.endswith('.gb') or infile.endswith('.gbk'):
            return is_genbank(infile)

        elif infile.endswith('.fa') or infile.endswith('.fasta') or infile.endswith('.faa') :
            return is_fasta(infile)

        else:
            print('Error: filetype of \'%s\' not recognized.' % (infile,))
            print('\tIf inputing a fasta file, please use .fa or .fasta')
            print('\tIf inputing a genbank file, please use .gb or .gbk')
            return False
    else:
        print('Error: file \'%s\'  does not exist' % (infile,))

def is_fasta(infile):
    g=SeqIO.parse(open(infile, 'r'), "fasta")
    res = peek(g)
    if res is None:
        return False
    else:
        return True

def is_genbank(infile):
    g=SeqIO.parse(open(infile, 'r'), "genbank")
    res = peek(g)
    if res is None:
        return False
    else:
        return True

def check_genbank_fmt(gbk):
    '''
    Takes a genbank file as input
    Checks if it has a locus of the proper length (in the future could
    check for other things?)
    If it does, returns the genbank.
    If not, rewrites locus name to the first 10 characters of the old one
    and returns an updated, correctly formated genbank.

    This is because the BioPython livrary does not allow genbanks with
    locus names greater than 10 characters. Another option is re-writing
    their function, because it honestly doesn't matter.
    '''
    count=0
    with open(gbk, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            count+=1
            if(len(str(record.name))) > 10:
                ofile = gbk.split('.')[0]+'.reformated.gbk'
                ofile_handle = open(ofile, 'w')
                print('\t\t##Warning##')
                print('Locus line on your gbk file is too long (common with ggKbase downloads)')
                print('\trewriting to first 10 characters of \'%s\'\n' % str(record.name))

                record.name = str(record.name)[:10]
                SeqIO.write(record, ofile_handle, 'genbank')
                ofile_handle.close()
            else:
                ofile = gbk

    if count >1:
        print('OH NO! multi-sequence genbank :( )')
        raise ValueError
    return ofile


def assign_nr_predictions(wd, gbk, nr_csv, pident, evalue, ofile):
    '''
    input = path to a gbk file
    output = path to a file with the NR additions


    Workflow:
    1) Fix the output CSV because sometimes there are commas in blast names
        (which messes up the pandas csv parser)
    2) Parse the CSV
    3) Choose best hit for each feature
    4) Separate into 'known' functional genes and 'unknowns'
    5) Add the functional gene anotations to an existing genbank file
    '''

    newcsv = nr_csv
    try:
        pd.read_csv(nr_csv)

    except pd.errors.ParserError:
        print('\tFound an error in your csv input, fixing...')
        newcsv = fix_psiblast(wd, nr_csv)
        print('\t\tFixed error')

    print('\tImporting csv into a pandas dataframe')
    print('\t...')
    df = parse_psiblast(newcsv, pident, evalue)

    print('\tChoosing the best hit for each feature')
    print('\t...')
    df1 = choose_hit(df)

    print('\tSorting out the hits to proteins with known functions')
    print('\t...')
    df2 = sort_unknowns(df1)
    feature_dict = return_feature_dict(df2[df2['funx']=='known'])

    print('\tAdding annotations to the genbank file')
    print('\t...')
    return add_NR_annotations(gbk, feature_dict, ofile)

def launch_psi_nr(faa, cluster= True):
    '''
    Update this so that it is passed a database and can run psiblast with any query/db
    '''
    cmdfile = 'nr_psi_search.sh'
    ocsv = faa.split('.')[0]+'.csv'
    if os.path.isfile(ocsv) and os.path.getsize(ocsv)>1:
        print('\tNR search found: \'%s\'' % (faaout,))
        print('\tTry running the \'assign\' workflow')

    else:
        scriptfile = open(cmdfile, 'w')
        if cluster:
            ofile = faa.split('.')[0]+'.gb'
            outfmt = '\'10 qseqid sseqid evalue bitscore length pident qcovs sblastname stitle scomname\''

            psi = 'psiblast -db /data1/blastdb/nr'
            psi = psi + ' -query '+faa
            psi = psi + ' -outfmt '+outfmt
            psi = psi + ' -num_threads 48 -out '+ocsv

            scriptfile.write(psi)

        scriptfile.close()
        cmd = 'qsub '+cmdfile
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

def launch_psi_envNR(faa, cluster=True):
    pass

def launch_psi_img(faa, cluster=True):
    pass

def find_psiblast_hmlgs(faa, cluster= True):
    '''
    Submits a job to the cluster searching the given faa file
    vs NR database

    Still needs some work, but is functional
    NOTE: ONLY works if you submit to the cluster. I do not know why
    '''
    launch_psi_nr(faa, cluster)
    launch_psi_envNR(faa, cluster)
    launch_psi_img(faa, cluster)

def run_fasta_workflow(**kwargs):
    if kwargs['faa'] != None:
        print('Protein predictions found.')
        faa = kwargs['faa']
    else:
        print('Checking for .faa file')
        faa = run_prodigal(kwargs['fasta'])

    print()
    print('Generating blank genbank file')
    print()
    gbk = generate_gbk(faa, kwargs['fasta'], kwargs['name'])

    print('Assinging functional proteins predicted by NR')
    print()
    assignedGbk = assign_nr_predictions(kwargs['wd'], gbk, kwargs['csv'],
                                        kwargs['pident'],kwargs['evalue'],kwargs['ofile'])
    print('DONE!')
    print('Functions assigned. See output genbank file: %s' % (assignedGbk,))

def run_genbank_workflow(**kwargs):
    checked_gbk = check_genbank_fmt(kwargs['gbk'])
    assignedGbk = assign_nr_predictions(kwargs['wd'], checked_gbk, kwargs['csv'],
                                        kwargs['pident'],kwargs['evalue'],kwargs['ofile'])
    print('DONE!')
    print('Functions assigned. See output genbank file: %s' % (assignedGbk,))

def argparser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #General Arguments
    GenArgs = parser.add_argument_group('GENERAL ARGUMENTS')
    GenArgs.add_argument('-h', action="help", help="show this help message and exit")
    GenArgs.add_argument('-i', "--input", help='path to input file', required=True)
    GenArgs.add_argument("-o", "--outdir", default = os.getcwd(),
            help='location of folder for output')
    GenArgs.add_argument("-a", "--proteins",
            help='path to protein predictions (must be done with prodigal)')
    GenArgs.add_argument("-d", "--dry_run", action = 'store_true',
            help='print commands that would be run [DOES NOT WORK]')
    GenArgs.add_argument("-psi", "--psioutput", required=True,
            help='path to outfile produced by psiblast search. Must be specifically formated.')
    GenArgs.add_argument('-f', "--filetype", help='Specify \'fasta\' or \'genbank\'')
    GenArgs.add_argument("-ogbk", "--outgbk",
            help='path final output file (a genbank)')


    Operational = parser.add_argument_group('OPTIONAL ARGUMENTS')
    Operational.add_argument('-p', '--pident',
            help='percent identity to define a homolog', default=30)
    Operational.add_argument('-e', '--evalue',
            help='evalue to define a homolog', default=0.0001)
    Operational.add_argument('-name', '--gbkname',
            help='Locus name for a genbank file. Must be <10 characters', default='123')

    args = parser.parse_args()
    return args

def main():
    '''
    Main controller
    '''
    a = argparser()
    pident = a.pident
    evalue = a.evalue

    faa = a.proteins
    csv = a.psioutput
    outdir = a.outdir
    infile = a.input
    name = a.gbkname
    wd = a.outdir
    ofile = a.outgbk

    if a.filetype != None:
        if a.filetype == 'fasta':
            print('Processing file as a fasta HERE')
            run_fasta_workflow(faa=faa,fasta=infile,csv=csv,wd=outdir,name=name,
                                pident=pident,evalue=evalue,ofile=ofile)

        if a.filetype == 'genbank':
            print('Processing file as a genbank')
            run_genbank_workflow(gbk=infile, wd=outdir,csv=csv,
                                pident=pident,evalue=evalue,ofile=ofile)

    else:
        if is_fasta(infile):
            print('Processing file as a fasta')
            run_fasta_workflow(faa=faa,fasta=infile,csv=csv,wd=outdir,name=name,
                                pident=pident,evalue=evalue,ofile=ofile)

        if is_genbank(infile):
            print('Processing file as a genbank')
            run_genbank_workflow(gbk=infile, wd=outdir,csv=csv,
                                pident=pident,evalue=evalue,ofile=ofile)


if __name__ == '__main__':
    main()
