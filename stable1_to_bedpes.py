from __future__ import print_function

import pandas as pd

'''
This script will parse an input gencode gtf and the supplementary table 1 from the INTEGRATE-Neo
paper, and will make a bedpe file for each sample in the table.

INPUTS:
    gencode.v22.annotation.gtf: Gencode v22 gtf annotation
    STable1.xlsx: Supplementary file from INTEGRATE-Neo

OUTPUTS:
    transgene_samples.txt: A list of all samples that had fusions (one per line)
    fusions/<TCGA ID>_fusions.bedpe for each <TCGA ID> in the table
'''

# A dict to store a mapping of genes to ENSEMBL IDs
gene_to_ensg = {}

with open('gencode.v22.annotation.gtf') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue    
        _line = line.strip().split('\t')
        if _line[2] == 'gene':
            __line = [x.strip() for x in _line[-1].split(';')]
            ensg_id = [x for x in __line if x.startswith('gene_id')][0].split()[1][1:-1]
            hugo_name = [x for x in __line if x.startswith('gene_name')][0].split()[1][1:-1]
            gene_to_ensg[hugo_name] = ensg_id.split('.')[0]

# Walk through each entry and write to an output file
outfile, off = '', None
infile = pd.read_excel('STable1.xlsx', header=0, skiprows=5)
samplefile = open('transgene_samples.txt', 'w')  # Stores sample names

try:
    for fusion in infile.itertuples():
        if fusion[1][:-3] not in outfile:
            # If a new samples is encountered, make a new outfile and close the old descriptor
            outfile = 'fusions/' + fusion[1][:-3] + '_fusions.bedpe'
            if off:
                off.close()
            off = open(outfile, 'w')
            print(fusion[1], file=samplefile)
        g1, g2= fusion[2].split('>>')
        print('chr%s\t' % fusion[3],
              '.\t',
              '%s\t' % fusion[4],
              'chr%s\t' % fusion[5],
              '%s\t' % fusion[6],
              '.\t',
              '%s-%s\t' % (gene_to_ensg[g1], gene_to_ensg[g2]),
              '.\t',
              '%s\t' % fusion[7], 
              '%s\t' % fusion[8],
              '.\t',
              '.\t',
              '%s\t' % g1,
              '%s' % g2, sep='', file=off)
finally:
    samplefile.close()
    off.close()