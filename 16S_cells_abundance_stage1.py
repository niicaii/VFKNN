
# -*- coding: utf-8 -*-

import os
import subprocess
import sys
import glob
import pandas as pd
import gzip
import sys
import re
# from .settings import logger
import argparse
# import gzip
# import threading
import warnings
warnings.filterwarnings('ignore')
# from .utils import buffer_count, simple_count
# from .settings import File, Setting, logger
# from .make_db import make_db

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
    help='Input file.')

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
    '-t',
    '--thread',
    metavar='INT',
    dest='thread',
    default=8,
    type=int,
    help='Number of threads. [8]')

required.add_argument(
    '--e1',
    metavar='FLOAT',
    dest='e1',
    default=1e-10,
    type=float,
    help='E-value cutoff for GreenGenes. [1e-10]')

required.add_argument(
    '--e2',
    metavar='FLOAT',
    dest='e2',
    default=3,
    type=float,
    help='E-value cutoff for Essential Single Copy Marker Genes. [3]')

required.add_argument(
    '--id',
    metavar='FLOAT',
    dest='id',
    default=45,
    type=float,
    help='Identity cutoff (in percentage) for Essential Single Copy Marker Genes. [45]')

required.add_argument(
    '--qcov',
    metavar='FLOAT',
    dest='qcov',
    default=0,
    type=float,
    help='Query cover cutoff (in percentage) for Essential Single Copy Marker Genes. [0]')

# required.add_argument(
#     '--database',
#     metavar='FILE',
#     dest='db',
#     default=None,
#     help='Customized database (indexed) other than SARG. [None]')

args=required.parse_args()

#%%
# green gene 85%
# ko.fasta, 19,951个氨基酸序列，都是些核糖体蛋白
# 共19,950行，给ko30.fasta的序列名分配了一个BA号

file = args.r_file 
# file = '/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/30soil/temp/pathogen/kraphlan/extractseq/10276_out.fa'
sample_name = str(args.sampname)
# sample_name = str(10276)
out = args.o_file
# out = '/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/30soil/temp/pathogen/kraphlan/'

ko30 = '/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/db/argoap/ko30'
ko30_structure = pd.read_table('/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/db/argoap/ko30_structure.txt', header=None, names=['sseqid', 'ko30'])
gg85 = '/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/db/argoap/gg85.fasta'
columns = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']

#%%

def buffer_count(file):
    '''
    Count number of lines of a fa/fq file, then divide by 2/4 to get reads.
    '''
    nlines = 0
    if re.search('.gz$', file):
        f = gzip.open(file, 'rt')
    else:
        f = open(file)

    char = f.read(1)  # get first character > or @
    f.seek(0)  # roll back

    buf_size = 1024 * 1024
    read_f = f.read

    buf = read_f(buf_size)
    while buf:
        nlines += buf.count('\n')
        buf = read_f(buf_size)

    f.close()
    if char == '>':
        return int(nlines / 2)
    elif char == '@':
        return int(nlines / 4)
    else:
        print('Unrecognized file format . Only fasta or fastq supported.')
        sys.exit(2)


def simple_count(file):
    '''
    Simple count bps of a file.
    '''
    nbps, nlines = 0, 0
    f = open(file)

    char = '>'
    record = False
    for line in f:
        if line[0] == char:
            nlines += 1
            record = True  # record next row

        elif record:
            nbps += len(line) - 1
            record = False

    f.close()
    return nbps, nlines

#%%
def uniq(df):
    '''
    uniq oap's blast best one
    '''
    df = df.sort_values(['qseqid', 'evalue', 'bitscore', 'length', 'pident'], ascending=[True, True, False, False, False])
    df = df.drop_duplicates(subset='qseqid', keep='first')
    
    # df1 = pd.DataFrame(columns=columns)._append(df[0:1])
    # num = int(len(df["qseqid"]))
    # for i in range(1, num):
    #     if df.loc[i, "qseqid"] == df1.loc[df1.index[-1], "qseqid"]:
    #         if df.loc[i, "bitscore"] <= df1.loc[df1.index[-1], "bitscore"]:
    #             if df.loc[i, "bitscore"] < df1.loc[df1.index[-1], "bitscore"]:
    #                 continue
    #             else:
    #                 if df.loc[i, "pident"] <= df1.loc[df1.index[-1], "pident"]:
    #                     if df.loc[i, "pident"] < df1.loc[df1.index[-1], "pident"]:
    #                         continue
    #                     else:
    #                         if df.loc[i, "evalue"] >= df1.loc[df1.index[-1], "evalue"]:
    #                             if df.loc[i, "evalue"] > df1.loc[df1.index[-1], "evalue"]:
    #                                 continue
    #                             else:
    #                                 if df.loc[i, "length"] <= df1.loc[df1.index[-1], "length"]:
    #                                     continue
    #                                 else:
    #                                     df1[-1:] = df[i:i + 1]
    #                         else:
    #                             df1[-1:] = df[i:i + 1]
    #                 else:
    #                     df1[-1:] = df[i:i + 1]
    #         else:
    #             df1[-1:] = df[i:i + 1]
    #     else:
    #         df1 = df1._append(df.loc[i])
    return df


#%%

def count_16s(file):
    '''
    Count 16S (GreenGenes 16S rRNA Database 85%) copy number using bwa (pre-filtering) and blastn (post-filtering).
    '''
    ## pre-filtering using bwa
    subprocess.run([
        'bwa', 'mem',
        '-t', str(args.thread),
        '-o', out + 'temp/tmp_16s_sam.16s.sam.tmp',
        gg85, file], check=True, stderr=subprocess.DEVNULL)

    tmp_16s_sam = out + 'temp/tmp_16s_sam.16s.sam.tmp'

    ## convert sam to fasta for later usage, note that reads can be duplicated
    with open(out + 'temp/tmp_16s_fa.16s.fa.tmp', 'w') as f:
        subprocess.run([
            'samtools',
            'fasta',
            '-F', '2308',
            tmp_16s_sam], check=True, stderr=subprocess.DEVNULL, stdout=f)

    tmp_16s_fa = out + 'temp/tmp_16s_fa.16s.fa.tmp'


    ## post-filter using blastn
    ## switch mt_mode if too little queries or too many threads, blast raises error if <2,500,000 bases per thread

    mt_mode = '1' if simple_count(tmp_16s_fa)[0] / args.thread >= 2500000 else '0'
    subprocess.run([
        'blastn',
        '-db', gg85,
        '-query', tmp_16s_fa,
        '-out', out + 'temp/tmp_16s_txt.16s.txt.tmp',
        '-outfmt', ' '.join(['6'] + columns),
        '-evalue', str(args.e1),
        '-max_hsps', '1',
        '-max_target_seqs', '1',
        '-mt_mode', mt_mode,
        '-num_threads', str(args.thread)], check=True, stderr=subprocess.DEVNULL)
    tmp_16s_txt = out + 'temp/tmp_16s_txt.16s.txt.tmp'


    ## process blastn results, store subject cover
    df = pd.read_table(tmp_16s_txt, header=None, names=columns)
    if len(df) == 0:
        print('No 16S-like sequences found in file.')
        # logger.warning(f'No 16S-like sequences found in file <{file.file}>.')
    else:
        if df['qseqid'].duplicated().sum() > 0:
            print('Duplicated sequences in 16S copy number calculation.')
            df = uniq(df)
            # subprocess.run([
            #     'python -m jcvi.formats.blast best -n 1 /media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/30soil/temp/pathogen/kraphlan/temp/tmp_16s_txt.16s.txt.tmp'
            #     ], 
            #     check=True, shell = True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            # print(result.stdout)
            # df = df[~df['qseqid'].duplicated()]
        return (df['length'] / df['slen']).sum()

#%%
def count_cells(file):
    '''
    Count Essential Single Copy Marker Genes (cell number) using diamond.
    '''
    ## filter using diamond
    subprocess.run([
        'diamond', 'blastx',
        '--db', f'{ko30}.dmnd',
        '--query', str(file),
        '--out', out + 'temp/tmp_cells_txt.cells.txt.tmp',
        '--outfmt', '6'] + columns + [
        '--evalue', str(args.e2),
        '--id', str(args.id),
        '--query-cover', str(args.qcov),
        '--max-hsps', '1',
        '--max-target-seqs', '1',
        '--threads', str(args.thread),
        '--quiet'], check=True, stderr=subprocess.DEVNULL)

    ## process blastx results, store subject coverage
    tmp_cells_txt = out + 'temp/tmp_cells_txt.cells.txt.tmp'

    df = pd.merge(pd.read_table(tmp_cells_txt, header=None, names=columns), ko30_structure, on='sseqid', how='left')

    if len(df) == 0:
        print('No marker-like sequences found in file.')
    else:
        if df['qseqid'].duplicated().sum() > 0:
            print('Duplicated sequences in cell number calculation.')
            df = uniq(df)
            # df = df[~df['qseqid'].duplicated()]

        return df.groupby('ko30').apply(lambda x: sum(x['length'] / x['slen'])).sum() / 30


# def extract_seqs(self, file):
#     '''
#     Prefilter ARGs, or other target sequences.
#     '''
#     ## diamond or bwa
#     if self.dbtype == 'prot':
#         subprocess.run([
#             'diamond', 'blastx',
#             '--db', f'{self.db}.dmnd',
#             '--query', file.file,
#             '--out', file.tmp_seqs_txt,
#             '--outfmt', '6', 'qseqid', 'full_qseq',
#             '--evalue', '10',
#             '--id', '60',
#             '--query-cover', '15',
#             '--max-hsps', '1',
#             '--max-target-seqs', '1',
#             '--threads', str(self.thread),
#             '--quiet'], check=True, stderr=subprocess.DEVNULL)
#     else:
#         subprocess.run([
#             'bwa', 'mem',
#             '-t', str(self.thread),
#             '-o', file.tmp_seqs_sam,
#             self.db, file.file], check=True, stderr=subprocess.DEVNULL)
#
#         with open(file.tmp_seqs_fa, 'w') as f:
#             subprocess.run([
#                 'samtools',
#                 'fasta',
#                 '-F', '2308',
#                 file.tmp_seqs_sam], check=True, stderr=subprocess.DEVNULL, stdout=f)
#
#         ## convert fa to tab to make results consistent
#         with open(file.tmp_seqs_txt, 'w') as w:
#             with open(file.tmp_seqs_fa) as f:
#                 for line in f:
#                     if line[0] == '>':
#                         qseqid = line[1:].rstrip().split(' ')[0]
#                     else:
#                         w.write(f'{qseqid}\t{line}')
#
#     ## give a new header for each target sequences, merge all sequences to a single file
#     with open(self.setting.extracted, 'a') as f:
#         df = pd.read_table(file.tmp_seqs_txt, header=None, names=['qseqid', 'full_qseq'])
#         if df['qseqid'].duplicated().sum() > 0:
#             logger.warning('Duplicated sequences in sequence extraction.')
#             df = df[~df['qseqid'].duplicated()]
#
#         cnt = 0
#         for qseqid, full_qseq in zip(df['qseqid'], df['full_qseq']):
#             f.write(f'>{file.sample_name}@{file.file_name}@{cnt}@{qseqid}\n{full_qseq}\n')
#             cnt += 1


#%%

def run(file):
    metadata = []
    metafile = sample_name + '_16S_cells_metadata.tsv'
    try:
        ## skip 16S/cells calculation if necessary
        nread = buffer_count(file)
        n16S = count_16s(file)
        ncell = count_cells(file)
        metadata.append([nread, n16S, ncell, sample_name])
    except subprocess.CalledProcessError:
        print(sample_name)
        print('Something is wrong with subprocess, skip.')

    pd.DataFrame(metadata, columns=['nRead', 'n16S', 'nCell', 'sample']).groupby('sample').sum().to_csv(out + '16S_cells_metadata/' + metafile, sep='\t')
    print('Finished.')

#%%

run(file)
