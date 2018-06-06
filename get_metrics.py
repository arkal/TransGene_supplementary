from __future__ import print_function
import dill
import gzip
import json
import os
import shutil
import subprocess

'''
This script will parse the outputs from `run_automated_transgene.py` assuming the outputs are
present as ./outputs . The results is a tsv matrix that produced supplementary table 1 from the
publication.

INPUTS:
        outputs/ : the outputs folder from `run_automated_transgene.py`
        sample_groups.dill: See `identify_sample_groups.py`
OUTPUTS:
        mutation_metrics.tsv: A TSV containing all results

'''

def split_call(call):
    pos = ''.join([x for x in call if x.isdigit()])
    return call[:call.find(pos)], pos, call[call.find(pos) + len(pos):]


with open('sample_groups.dill') as df:
    samples = dill.load(df)

with open('mutation_metrics.tsv', 'w') as off:
    print('Sample',
          'Total_input_muts',
          'Total_PASS_muts',
          'Total_accepted_muts',
          'Total_PASS_snvs',
          'Total_accepted_snvs',
          'Total_PASS_indels',
          'Total_accepted_indels',
          'Total_fusions',
          'Accepted_fusions',
          '9_mer_snvs',
          '9_mer_indels',
          '9_mer_fusions',
          '10_mer_snvs',
          '10_mer_indels',
          '10_mer_fusions',
          '15_mer_snvs',
          '15_mer_indels',
          '15_mer_fusions',
          'time_taken',
          sep='\t', file=off)

    for sample in samples:
        total_muts = pass_muts = accepted_muts = _total_muts = 0
        total_fusions = accepted_fusions = 0
        peptides = {9: {}, 10: {}, 15: {}}
        # Read the statistics produced by run_automated_transgene.py
        with open('outputs/%s/%s.log2' % (sample, sample)) as lf:
            total_muts = int(lf.readline().split(':')[1].strip())
            pass_muts = int(lf.readline().split(':')[1].strip())
            time_taken = round(float(lf.readline().split(':')[1].strip()), 2)
        # Parse the transgened vcf
        input_pass_snvs = input_pass_indels = input_accepted_snvs = input_accepted_indels = 0
        with open('outputs/%s/%s_transgened.vcf' % (sample, sample)) as ivf:
            for line in ivf:
                if not line.startswith('#'):
                    _total_muts += 1
                    line = line.strip().split()
                    if len(line[3]) == len(line[4]):
                        input_pass_snvs += 1
                    else:
                        input_pass_indels += 1                    
                    if line[6] == 'PASS':
                        accepted_muts += 1
                        if len(line[3]) == len(line[4]):
                            input_accepted_snvs += 1
                        else:
                            input_accepted_indels += 1

        assert pass_muts == _total_muts
        assert input_pass_snvs + input_pass_indels == _total_muts
        assert input_accepted_snvs + input_accepted_indels == accepted_muts
        # Parse the input fusions
        with open('fusions/%s_fusions.bedpe' % sample) as iff:
            for line in iff:
                total_fusions +=1
        # Parse the accepted fusions
        with open('outputs/%s/%s_transgened.bedpe' % (sample, sample)) as iff:
            for line in iff:
                accepted_fusions +=1
        # count peptides
        for n in [9, 10, 15]:
            if os.path.exists('outputs/%s/%s_tumor_%s_mer_snpeffed.faa.map' % (sample, sample, n)):
                ivf = open('outputs/%s/%s_tumor_%s_mer_snpeffed.faa.map' % (sample, sample, n))
            elif os.path.exists('outputs/%s/%s_tumor_%s_mer_peptides.faa.map' % (sample, sample, n)):
                ivf = open('outputs/%s/%s_tumor_%s_mer_peptides.faa.map' % (sample, sample, n))
            else:            
                peptides[n]['f'] = 0
                peptides[n]['s'] = 0
                peptides[n]['i'] = 0
                continue
            temp = json.load(ivf)
            ivf.close()
            peptides[n]['f'] = len([x_ for x_ in temp if 'FUSION' in temp[x_]])
            peptides[n]['i'] = 0
            non_fusions = len(temp) - peptides[n]['f']

            for k, v in temp.items():
                if 'FUSION' in v:
                    continue
                v = v.split('\t')
                v_ = v[-1].split(',')
                for v__ in v_:
                    v__ = split_call(str(v__).split('#')[0].split('_')[1])
                    if '?' in v__[2] or len(v__[0]) != len(v__[2]):
                        peptides[n]['i'] += 1
                        break
            peptides[n]['s'] = non_fusions - peptides[n]['i']
        print(sample,
              total_muts,
              pass_muts,
              accepted_muts,
              input_pass_snvs,
              input_accepted_snvs,
              input_pass_indels,
              input_accepted_indels,
              total_fusions,
              accepted_fusions,
              peptides[9]['s'],
              peptides[9]['i'],
              peptides[9]['f'],
              peptides[10]['s'],
              peptides[10]['i'],
              peptides[10]['f'],
              peptides[15]['s'],
              peptides[15]['i'],
              peptides[15]['f'],
              time_taken,
              sep='\t', file=off)
        
