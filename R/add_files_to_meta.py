#!/usr/bin/env python
# coding: utf-8
"""
Add file location to metadata 
"""

from pathlib import Path
import pandas as pd
import numpy as np
from itertools import chain

df = pd.read_excel(
    '/prj/Niels_Gehring/nmd_transcriptome/nmd_transcriptome_metadata.xlsx',
    engine='openpyxl')
df = df.iloc[:62]

p1 = Path('/prj/Niels_Gehring/batchMar_19W14/workflow/mapping')
p2 = Path('/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batchMar2_19W14/workflow/mapping')
p3 = Path('/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batch3_19W14/workflow/mapping')
p4 = Path('/prj/Niels_Gehring/DFG_seq_Nanopore/bams')
ext = '*L002_STARmapping/Aligned.noS.bam'


bams = list(
    chain(
        p1.glob(ext), 
        p2.glob(ext),
        p3.glob(ext),
        p4.glob("*_DRS_all_pass.bam")
    )
)

bams = np.array(list(bams)).astype(str)


def find_match(a, b):
    match = np.char.find(a, b)
    match = np.flatnonzero(match != -1)
    if match.shape[0]==1:
        return match[0]
    else:
        return None

matches = [find_match(bams, x) for x in df['CCG_Sample_ID']]
# sanity check
pd.isnull(matches).any()
df.loc[pd.isnull(matches).nonzero()[0]]
df.loc[:, 'bams'] = pd.Series(bams).reindex(matches).values
# none should be missing
df['bams'].apply(lambda x: Path(x).exists()).all()

def find_files(paths, ext):
    for p in paths:
        yield from Path(p).glob(ext)


p1 = Path('/prj/Niels_Gehring/batchMar_19W14/workflow/raw_reads')
p2 = Path('/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batchMar2_19W14/workflow/raw_reads')
p3 = Path('/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batch3_19W14/workflow/raw_reads')
ext = '*_R1_001.fastq.gz'


def find_files(paths, ext):
    for p in paths:
        yield from Path(p).glob(ext)

ext = '*_R1_001.fastq.gz'
paths = [
    "/prj/Niels_Gehring/batchMar_19W14/workflow/raw_reads",
    "/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batchMar2_19W14/",
    "/prj/Niels_Gehring/KD_different_cell_lines_with_SMG6_SMG7/batch3_19W14/"    
]

raw_reads = find_files(paths, ext)
raw_reads = [x for x in raw_reads if x.exists()]
raw_reads = [str(x) for x in raw_reads]
raw_reads_match = [find_match(raw_reads, x) for x in df['CCG_Sample_ID']]
df['raw_reads'] = pd.Series(raw_reads).reindex(raw_reads_match).values
df.to_csv("/prj/Niels_Gehring/nmd_transcriptome/phase2/metadata_w_files.csv")

