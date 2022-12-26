import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EXPANSION_SIZE = 350

def repeat_to_length(s, wanted):
    '''
    Source: https://stackoverflow.com/questions/3391076/repeat-string-to-certain-length
    '''
    return (s * (wanted//len(s) + 1))[:wanted]

# Load the chr1 reference genome sequence
chr1 = SeqIO.read(snakemake.input['ref'], "fasta")
chr1_new = SeqRecord(Seq(''), name='chr1', id='chr1')

repeats = pd.read_csv(snakemake.input['repeats'])

for i, repeat in repeats.iterrows():
    start = repeat['start']
    stop = repeat['stop']
    motif = repeat['motif']

    if i > 0:
        prev_stop = repeats.loc[i - 1, 'stop']
    else:
        prev_stop = 0
    
    if i < len(repeats) - 1:
        next_start = repeats.loc[i + 1, 'start']
    else:
        next_start = len(chr1)

    # Q: Do we need to lower or upper here?
    seq_to_insert = repeat_to_length(motif, EXPANSION_SIZE).lower()

    # Modify the chr1 sequence by expanding the repeat
    chr1_new.seq += chr1.seq[prev_stop:start] + seq_to_insert + chr1.seq[stop:next_start]

# Create a reference for the expansion and flanking region around it
# record = SeqRecord(chr1_new)

print('CHR LENGTH:', len(chr1_new.seq))

SeqIO.write(chr1_new, snakemake.output[0], 'fasta')
