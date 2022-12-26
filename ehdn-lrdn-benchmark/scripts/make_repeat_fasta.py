import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EXPANSION_SIZE = 350
FLANK_SIZE = 1000

def repeat_to_length(s, wanted):
    '''
    Source: https://stackoverflow.com/questions/3391076/repeat-string-to-certain-length
    '''
    return (s * (wanted//len(s) + 1))[:wanted]

# Load the chr1 reference genome sequence
chr1 = SeqIO.read(snakemake.input['ref'], "fasta")

repeats = pd.read_csv(snakemake.input['repeats']).set_index('index')
start = repeats.loc[snakemake.wildcards['repeat'], 'start']
stop = repeats.loc[snakemake.wildcards['repeat'], 'stop']
motif = repeats.loc[snakemake.wildcards['repeat'], 'motif']

# Q: Do we need to lower or upper here?
seq_to_insert = repeat_to_length(motif, EXPANSION_SIZE).lower()

# Modify the chr1 sequence by expanding the repeat
chr1.seq = chr1.seq[:start] + seq_to_insert + chr1.seq[stop:]

stop = start + EXPANSION_SIZE

# Create a reference for the expansion and flanking region around it
record = SeqRecord(chr1.seq[start - FLANK_SIZE:stop + FLANK_SIZE])

SeqIO.write(record, snakemake.output[0], 'fasta')
