# -*- coding: utf-8 -*-

import os
import subprocess
import sys
import pandas as pd
import argparse
import warnings
warnings.filterwarnings('ignore')

#%%
"""
添加参数模块
"""

required = argparse.ArgumentParser(description='manual to this script')

required.add_argument(
    '-i',
    '--input',
    metavar='FILE',
    dest='r_file',
    required=True,
    help='Input blastfile,such as vf blastout.')

required.add_argument(
    '-m',
    '--meta',
    metavar='FILE',
    dest='meta',
    required=True,
    help='Input metafile preduced in stage 1.')

required.add_argument(
    '-o',
    '--out',
    dest='o_file',
    required=True,
    help='Output dir.')

required.add_argument(
    '--sampname',
    dest='sampname',
    required=True,
    help='which sample and its name.')

required.add_argument(
    '--dbtype',
    dest='dbtype',
    required=True,
    help='which db type ,prot or nucl.')

required.add_argument(
    '-t',
    '--thread',
    metavar='INT',
    dest='thread',
    default=8,
    type=int,
    help='Number of threads. [8]')

required.add_argument(
    '--e',
    metavar='FLOAT',
    dest='e',
    default=1e-5,
    type=float,
    help='E-value cutoff for target sequence. [1e-5]')

required.add_argument(
    '--id',
    metavar='FLOAT',
    dest='id',
    default=80,
    type=float,
    help='Identity cutoff (in percentage) for target sequence. [80]')

required.add_argument(
    '--qcov',
    metavar='FLOAT',
    dest='qcov',
    default=75,
    type=float,
    help='Query cover cutoff (in percentage) for  target sequence. [75]')

required.add_argument(
    '--length',
    metavar='INT',
    default=25,
    type=int,
    help='Aligned length cutoff (in amino acid) for target sequences. [25]')

# required.add_argument(
#     '--database',
#     metavar='FILE',
#     dest='db',
#     default=None,
#     help='Customized database (indexed) other than SARG. [None]')

args=required.parse_args()
#%%
blastout = str(args.r_file)
# blastout = 'E:/bioinformatics_linux/30soil/temp/pathogen/kraphlan/etractseqblastvf/10276_vfblastout.txt'
sample_name = str(args.sampname)
# sample_name = str(10276)
out = args.o_file
# out = 'E:/bioinformatics_linux/30soil/temp/pathogen/kraphlan/temp/'
metadata = str(args.meta)
# metadata = 'E:/bioinformatics_linux/30soil/temp/pathogen/kraphlan/16S_cells_metadata/10280_16S_cells_metadata.tsv'
dbtype = str(args.dbtype)

columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'bitscore']

# columnsout = ['sseqid', 'qseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore', 'qcov', 'scov']

#%%

def __init__(self, options):

    ## if zero entries in metadata then cannot normalize
    metadata = pd.read_table(str(args.meta), dtype={'sample': str})
    if (metadata[['nRead', 'n16S', 'nCell']] == 0).any(axis=None):
        print('Found zero reads/16s/cells in metadata file. These samples will be ignored.')
        metadata = metadata[~(metadata == 0).any(axis=1)]

    if len(metadata) == 0:
        print('No valid sample in metadata file. No further normalization will be made.')
        sys.exit(2)


#%%
def tpm():
    '''
    Join extracted target sequences with metadata and structures. Aggregate according to levels (type/subtype/gene).
    '''
    # print('Merging files ...')
    df = pd.read_table(blastout, header=None, names=columns, dtype={'sseqid': str})

    if len(df) == 0:
        print('No target sequence detected in file , no further normalization will be made.')
        sys.exit(2)

    if df.isnull().any(axis=None):
        print('blast file cannot be parsed. Please check BLAST output file (--blastout).')
        sys.exit(2)

    ## further filtering
    # 此处
    qcov = args.qcov / 100 / 3 if dbtype == 'prot' else args.qcov / 100  # for prot need to lower qcov as qlen is in bp # qcov 75 指的核酸比核酸
    length = args.length if dbtype == 'prot' else args.length * 3  # for nucl need to convert aa cut to bp # 25 Aligned length cutoff (in amino acid 氨基酸) for target sequences.

    df['qcov'] = df['length'] / df['qlen']
    df['scov'] = df['length'] / df['slen']

    ## filter
    df = df[(df['pident'] >= args.id) & (df['evalue'] <= args.e) & (df['length'] >= length) & (df['qcov'] >= qcov)]

    ## deduplicate
    df = df.sort_values(['qseqid', 'evalue', 'bitscore', 'length', 'pident'], ascending=[True, True, False, False, False])
    df = df.drop_duplicates(subset='qseqid', keep='first')

    ## tpm
    df['nreads_within_s'] = df.groupby('sseqid')['qseqid'].transform('count')
    df['averl_within_s'] = df.groupby('sseqid')['length'].transform('mean')
    df['rpk_within_s'] = (df['nreads_within_s']*df['averl_within_s']) / (df['slen']/1000)
    df['tpm_within_s'] = (df['rpk_within_s']*1e6)/sum(df.groupby('sseqid')['rpk_within_s'].transform('mean'))
    # df.to_csv(out + 'temp/vftpmtemp/' + str(args.sampname) + '_vfblast_filter.tsv', sep='\t')
    # dft = df.groupby('sseqid')['tpm'].agg('mean').to_frame()
    meta = pd.read_table(metadata, header=0, dtype={'sample': str})
    dft = df.assign(taxid=meta.loc[0,'sample'],
              nReads=meta.loc[0,'nRead'],
              n16S=meta.loc[0,'n16S'],
              nCell=meta.loc[0,'nCell'])
    return dft
#%%
def run():
    df = tpm()
    outfile = 'vftpm/' + str(args.sampname) + '_vf_tpm.tsv'
    df.to_csv(out + outfile, sep='\t', header=True, index=False)
    print(str(args.sampname))
    print('Finished.')

#%%
run()