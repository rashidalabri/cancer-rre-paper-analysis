import json
import random

import pandas as pd

random.seed(1)

N_EXPANSIONS = 10

# Load variants file
with open(snakemake.input[0]) as f:
    variants = json.load(f)

# Extract simple repeats in chr1
variants = [variant for variant in variants if isinstance(
    variant['ReferenceRegion'], str) and variant['ReferenceRegion'].startswith('chr1:')]

# Sample a set number of repeats
variants = random.sample(variants, k=N_EXPANSIONS)

df = {'index': [], 'chr': [], 'start': [], 'stop': [], 'motif': []}
for v in variants:
    chr, loc = v['ReferenceRegion'].split(':')
    start, stop = loc.split('-')
    index = f'{chr}_{start}_{stop}'
    motif = v['LocusStructure'][1:-2]

    df['index'].append(index)
    df['chr'].append(chr)
    df['start'].append(start)
    df['stop'].append(stop)
    df['motif'].append(motif)

df = pd.DataFrame.from_dict(df)

df.to_csv(snakemake.output[0], index=False)
