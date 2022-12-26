import pandas as pd
from Bio.Seq import Seq
import subprocess
import json


SAMTOOLS_PATH = '/scg/apps/software/samtools/1.11/bin/samtools'


window = 1000

profiles = [pd.read_table(f'../data/output/{i}.locus.tsv') for i in range(0, 11)]

repeats = '../data/input/repeats.csv'
repeats = pd.read_csv(repeats)

ranges = {(row['start'], row['stop']): row['motif'] for _, row in repeats.iterrows()}

# Make motif equivalence set
def generate_set(motif):
    n = len(motif)
    result = set()
    for i in range(n):
        seq = Seq(motif[i:]+motif[:i])
        result.add(str(seq))
        result.add(str(seq.reverse_complement()))
#         result.add(str(seq) + str(seq))
#         result.add(str(seq.reverse_complement) + str(seq.reverse_complement))
    return result

equiv_set = {}
for _, motif in ranges.items():
    equiv_set[motif] = generate_set(motif)
    
# print(equiv_set['TAAA'])
# exit()

# extracted = pd.DataFrame(columns=profile.columns)

def get_total_read_depth(bam_file_path, chrom, start, stop):
    region = '{}:{}-{}'.format(chrom, start, stop)

    samtools_cmd = '{} depth -r {} {}'.format(SAMTOOLS_PATH, region, str(bam_file_path))
    awk_cmd = ['awk', '{d+=$3}END{print d}']

#     print(samtools_cmd, flush=True)

    ps = subprocess.Popen(samtools_cmd.split(' '), stdout=subprocess.PIPE)
    result = subprocess.check_output(awk_cmd, stdin=ps.stdout)
    ps.wait()

    return int(result.strip())

result_df = pd.DataFrame(columns=['position', 'motif', 'amp', 'glob_irr', 'local_irr'])

for i, p in enumerate(profiles):
    print(f'--- [copy no. {i}] ---')
    
    f = open(f'../data/output/{i}.str_profile.json')
    data = json.load(f)
    global_read_depth = data['Depth']
    f.close()
        
    found = False
    for _, row in p.iterrows():
        for (start, end), motif in ranges.items():
            if start - window < row['start'] and row['start'] < start + window and end - window < row['end'] and row['end'] < end + window:
                if str(row['motif']) in equiv_set[motif]:
                    found = True
                    irr = row['num_anc_irrs']
                    glob_irr = row['norm_num_anc_irrs']
                    
                    buff_start = int(row['start']) - 500
                    buff_stop = int(row['end']) + 500
                    local_read_count = get_total_read_depth(f'../data/input/sorted/{i}.bam', row['contig'], buff_start, buff_stop)
                    local_read_depth = local_read_count / (buff_stop - buff_start)
                    local_irr = irr / local_read_depth
                    
                    result_df = result_df.append({
                        'position': f'chr1:{start}-{end}',
                        'motif': motif,
                        'amp': i,
                        'glob_avg_read_depth': global_read_depth,
                        'local_avg_read_depth': local_read_depth,
                        'glob_irr': glob_irr,
                        'local_irr': local_irr,
                    }, ignore_index=True)
                    
                    
                    
                    print(f'Found {motif};{start}-{end}, Global Norm. IRR: {glob_irr}, Local Norm. IRR: {local_irr}')
    if not found:
        print('No expansions recovered.')
               

result_df.to_csv('../data/output/results.tsv', index=False, sep='\t')