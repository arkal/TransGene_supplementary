from __future__ import print_function

import dill


'''
This script will parse an input gdc sample sheet for vcfs, the transgene_samples.txt file from 
`stable1_to_bedpes.py` (samples with fusions) to produce a dill file that contains the gdc UUIDs
for the RNA bam and the VCF for each sample.

INPUTS:
    gdc_sample_sheet_vcfs.tsv: Sample sheet for DNA bams
    transgene_samples.txt: See `stable1_to_bedpes.py`

OUTPUTS:
    sample_groups.dill: A dill dump of a dictionary with structure

                {TCGA_ID: {rna: GDC_UUID, vcf: GDC_UUID}}
'''

def read_gdc_manifest(manifest):
    output = {}
    with open(manifest) as f:
        for line in f:
            if line.startswith('File'):
                continue
            else:
                line = line.strip().split('\t')
                output[line[5].split(',')[0]] = line[1].split('.')[0]
    return output

def get_shortlisted_samples(filename):
    output = []
    with open(filename) as f:
        for line in f:
            output.append(line.strip())
    return output


manifest = read_gdc_manifest('gdc_sample_sheet_vcfs.tsv')
samples = get_shortlisted_samples('transgene_samples.txt')

output = {x: {'vcf': manifest[x]} for x in samples}

with open('/Users/arjun/paper_work/protect/input_files/gdc_sample_sheet_rna_bams_filtered.tsv') as iff:
    for line in iff:
        _line = line.strip().split('\t')
        if _line[5] in output:
            output[_line[5]]['rna'] = _line[0]

output = {x: y for x, y in output.items() if len(y) == 2}

with open('sample_groups.dill', 'w') as df:
    dill.dump(output, df)