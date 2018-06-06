from __future__ import print_function
import dill
import gzip
import json
import os
import re
import shutil
import subprocess


with open('sample_groups.dill') as df:
    samples = dill.load(df)

regex = re.compile('#[^,]+#')
for n in 9,10,15:
    print(n)
    for sample in samples:
      # count peptides
      if os.path.exists('outputs/%s/%s_tumor_%s_mer_snpeffed.faa.map' % (sample, sample, n)):
          with open('outputs/%s/%s_tumor_%s_mer_snpeffed.faa.map' % (sample, sample, n)) as ivf:
              temp = json.load(ivf)
              if any(re.findall(regex, y)!=[] for y in temp.values()):
                print(sample)
      elif os.path.exists('outputs/%s/%s_tumor_%s_mer_peptides.faa.map' % (sample, sample, n)):
          with open('outputs/%s/%s_tumor_%s_mer_peptides.faa.map' % (sample, sample, n)) as ivf:
              temp = json.load(ivf)
              if any(re.findall(regex, y)!=[] for y in temp.values()):
                print(sample)
      else:            
          raise


# tree -fi | grep log$ | while read line; do echo $line; grep "Merging" ${line}; done
# will give number of times the mutants affected the same codon
