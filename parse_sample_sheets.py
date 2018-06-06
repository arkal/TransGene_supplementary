from __future__ import print_function

import dill

from collections import defaultdict

'''
This script will parse input gdc samples sheets for DNA and RNA-Seq bam files to only retain files
corresponding to samples in an input protect_samples.txt file. The three files are expected to be in
the same directory as this script. The files retained will futher be filtered to only retain one 
DNA tumor bam, DNA normal bam, and one RNA tumor bam for each sample.

INPUTS:
    protect_samples.txt: a text file with one TCGA ID per line
    gdc_sample_sheet_dna_bams.tsv: Sample sheet for DNA bams
    gdc_sample_sheet_rna_bams.tsv: Sample sheet for RNA bams

OUTPUTS:
    gdc_sample_sheet_dna_bams_filtered.tsv: Filtered sample sheet for DNA bams
    gdc_sample_sheet_rna_bams_filtered.tsv: Filtered sample sheet for RNA bams
'''

samples = set()
with open('protect_samples.txt') as sf:
    for line in sf:
        samples.add(line.strip())

# Stores all files associated with each seq type for each sample.
sample_files = {x: {'dna': [], 'rna': []} for x in samples}

for nuc in 'dna', 'rna':
    with open('gdc_sample_sheet_%s_bams.tsv' % nuc) as iff:
        for line in iff:
            _line = line.strip().split('\t')
            if _line[5] in samples:
                sample_files[line.strip().split('\t')[5]][nuc].append((_line[6], _line[0]))

# A dict of all sample groups that contain the required files to complete proper trios.
out_groups = {x:y for x,y in sample_files.items() if len(y['rna'])==1 and len(y['dna'])==2}

# A dict of all sample group that contains more or less than the required files to complete proper trios.
fix_groups = {x:y for x,y in sample_files.items() if len(y['rna'])!=1 or len(y['dna'])!=2}

# For groups that have too many RNA files, discard the normal bam files.
for sample in fix_groups:
    if len(fix_groups[sample]['rna']) != 1:
        fix_groups[sample]['rna'] = [x for x in fix_groups[sample]['rna'] if x[0].endswith(('01A', '01B'))]

# Another dict of all sample groups that contain the required files to complete proper trios.
fixed_groups_1 = {x:y for x,y in fix_groups.items() if len(y['rna'])==1 and len(y['dna'])==2}

# Update to remove fixed groups
fix_groups = {x:y for x,y in fix_groups.items() if len(y['rna'])!=1 or len(y['dna'])!=2}

# make sure that the tumor dna and rna come from the same source (01/-01B)
for sample in fix_groups:
    if len(fix_groups[sample]['dna']) != 2:
        assert len(fix_groups[sample]['dna']) > 2
        assert len({x[0].split('-')[-1] for x in fix_groups[sample]['dna']}) == len(fix_groups[sample]['dna'])
        outfiles = {x[0].split('-')[-1]: x for x in fix_groups[sample]['dna']}
        fix_groups[sample]['dna'] = [outfiles['01B'] if '01B' in outfiles else outfiles['01A']]
        fix_groups[sample]['dna'].append(outfiles['11A'] if '11A' in outfiles else outfiles['10A'])
        assert len(fix_groups[sample]['dna']) == 2

# Another dict of all sample groups that contain the required files to complete proper trios.
fixed_groups_2 = {x:y for x,y in fix_groups.items() if len(y['rna'])==1 and len(y['dna'])==2}

# Update to remove fixed groups
fix_groups = {x:y for x,y in fix_groups.items() if len(y['rna'])!=1 or len(y['dna'])!=2}

assert not fix_groups, 'There are unprocessed samples.'

out_groups.update(fixed_groups_1)
out_groups.update(fixed_groups_2)

# Sanity checks
for group in out_groups:
    assert out_groups[group]['rna'][0][0] in (out_groups[group]['dna'][0][0], out_groups[group]['dna'][1][0])
    assert {x[0].split('-')[-1][:-1] for x in out_groups[group]['dna']} in [{'01', '10'}, {'01', '11'}]

# Print results
for nuc in 'dna', 'rna':
    group_samples = set()
    for sample in out_groups:
        group_samples.update({x[1] for x in out_groups[sample][nuc]})
    with open('gdc_sample_sheet_%s_bams.tsv' % nuc) as iff, open('gdc_sample_sheet_%s_bams_filtered.tsv' % nuc, 'w') as off:
        print(iff.readline(), end='', file=off)
        for line in iff:
            _line = line.strip().split('\t')
            if _line[0] in group_samples:
                print(line, end='', file=off)
