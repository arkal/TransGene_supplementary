#!/usr/bin/env python2.7
from __future__ import print_function
from toil.job import Job

import argparse
import dill
import gzip
import os
import shutil
import subprocess
import time

'''
This program will run transgene on every sample in the input `sample_groups.dill` dictionary. The
RNA bama nd the vcf are downloaded for each sample and then deleted at the end of the run to prevent
unnecessary storage of the temporary files.

To run, try 
        
                    python run_automated_transgene.py --help
'''

def run_transgene(job, sample, sample_info, common, creds, output_folder):
    """
    Get the files from GDC using the credentials file. Run Transgene on it.
    """
    total_muts = pass_muts = 0
    # Get the vcf
    parameters = ['/usr/local/bin/gdc-client',
                  'download',
                  '-n', '15',
                  '-t', creds,
                  sample_info['vcf']]
    subprocess.check_call(parameters)
    in_vcf_file = [f for f in os.listdir(sample_info['vcf']) if f.endswith('vcf.gz')][0]
    in_vcf_file = os.path.join(os.getcwd(), sample_info['vcf'], in_vcf_file)
    # Parse the vcf
    with gzip.open(in_vcf_file) as ivf, open(sample + '.vcf', 'w') as ovf:
        for line in ivf:
            if line.startswith('#'):
                print(line.strip(), file=ovf)
            else:
                total_muts += 1
                if line.strip().split()[6] == 'PASS':
                    pass_muts += 1
                    print(line.strip(), file=ovf)
    vcf_file = os.path.join(os.getcwd(), ovf.name)
    # Get the rna bam
    parameters = ['/usr/local/bin/gdc-client',
                  'download',
                  '-t', creds,
                  sample_info['rna']]
    subprocess.check_call(parameters)
    rna_bam = [f for f in os.listdir(sample_info['rna']) if f.endswith('bam')][0]
    rna_bam = os.path.join(os.getcwd(), sample_info['rna'], rna_bam)
    # rename the bai file for now
    os.rename(rna_bam[:-1] + 'i', rna_bam + '.bai')
    # Get the bedpe
    fusions = os.path.join(common['fusion_folder'], sample + '_fusions.bedpe')

    # Create the output directory
    os.mkdir(os.path.join(output_folder, sample))

    # Run Transgene
    parameters = ['transgene',
                  '--prefix', os.path.join(output_folder, sample, sample),
                  # SNV + INDEL
                  '--vep', vcf_file, 
                  '--indel_extend_length', '30',
                  # FUSION
                  '--fusions', fusions,
                  '--filter_mt_fusions',
                  '--filter_ig_pairs',
                  '--filter_rna_gene_fusions',
                  '--filter_readthroughs',
                  # RNA
                  '--rna_file', rna_bam,
                  # AUX
                  '--genome', common['genome'],
                  '--peptides', common['peptides'], 
                  '--annotation', common['annotation'],
                  '--transcripts', common['transcripts'],
                  # MISC
                  '--log_file', os.path.join(output_folder, sample, sample + '.log'), 
                  '--cores', '15'
    ]
    start = time.time()
    subprocess.check_call(parameters)
    time_taken = time.time() - start
    with open(os.path.join(output_folder, sample, sample + '.log2'), 'w') as lf:
        print('Total_muts: %s\nPASS_muts: %s\ntime_taken:%s' % (total_muts, pass_muts, time_taken), file=lf)


def launchpad(job, sample_groups_file, input_folder, output_folder, creds):
    """
    The start of the path.
    """
    # Load the sample groups
    with open(sample_groups_file, 'r') as df:
        sample_groups = dill.load(df)

    # Make the common dict
    common = {'fusion_folder': os.path.join(input_folder, 'fusions'),
              'peptides': os.path.join(input_folder, 'gencode.v22.pc_translations.fa'),
              'transcripts': os.path.join(input_folder, 'gencode.v22.pc_transcripts.fa'),
              'annotation': os.path.join(input_folder, 'gencode.v22.annotation.gtf'),
              'genome': os.path.join(input_folder, 'GRCh38.d1.vd1.fa')
    }

    for sample, sample_info in sample_groups.items():
        job.addChildJobFn(run_transgene,
                          sample,
                          sample_info,
                          common,
                          creds,
                          output_folder,
                          cores=15)


def main():
    """
    This is the main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_groups', '-S', dest='sample_groups', help='sample_groups.dill', 
                        type=str, required=True)
    parser.add_argument('--creds', '-C', dest='creds', help='GDC token file.', type=str,
                        required=True)
    parser.add_argument('--output_folder', '-O', dest='output_folder', help='Output folder.',
                        type=str, required=False, default='outputs')
    parser.add_argument('--input_folder', '-I', dest='input_folder', help='Input folder.',
                        type=str, required=False, default=os.getcwd())
    Job.Runner.addToilOptions(parser)
    params = parser.parse_args()

    params.sample_groups = os.path.abspath(params.sample_groups)
    params.creds = os.path.abspath(params.creds)
    params.output_folder = os.path.abspath(params.output_folder)
    params.input_folder = os.path.abspath(params.input_folder)

    start = Job.wrapJobFn(launchpad, params.sample_groups, params.input_folder,
                          params.output_folder, params.creds, cores=1)
    Job.Runner.startToil(start, params)
    return None

if __name__ == '__main__':
    main()