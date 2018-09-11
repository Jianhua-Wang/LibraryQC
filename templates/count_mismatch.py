#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import re

if len(sys.argv) == 1:
    sys.argv=["count_mismatch.py","$sam"]

sam = sys.argv[1]
mismatch = 'mismatch.txt'
pbc_qc_txt = 'pbc_qc.txt'

with open(mismatch, 'w') as f:
    f.write('read_id\\tlib_id\\tsample_seq\\tmis\\tMAPQ\\n')
    with open(sam,'r') as fp:
        ith = 0
        for line in fp:
            if  re.match(r'^@',line):
                continue
            ith += 1
            line_sp = line.split()
            if line_sp[1] == '0':
                f.write(line_sp[0]+'\\t'+line_sp[2]+'\\t'+line_sp[9]+'\\t'+line_sp[12][-1]+'\\t'+line_sp[4]+'\\n')
            else:
                f.write(line_sp[0]+'\\t*\\t'+line_sp[9]+'\\tunmapped\\t'+line_sp[4]+'\\n')

total_df = pd.read_table(mismatch)[['lib_id','sample_seq','mis']].drop_duplicates(subset = 'sample_seq')
distinct_df = total_df.groupby('lib_id').size()
value_counts = distinct_df.value_counts()

pbc_qc = pd.Series()
pbc_qc['Item'] = 'Value'
pbc_qc['TotalReadPairs'] = len(total_df)
pbc_qc['DistinctReadPairs'] = pbc_qc['TotalReadPairs'] - distinct_df['*']
pbc_qc['OneReadPair'] = value_counts[1]
pbc_qc['TwoReadPairs'] = value_counts[2]
pbc_qc['NRF'] = pbc_qc['DistinctReadPairs']/pbc_qc['TotalReadPairs']
pbc_qc['PBC1'] = pbc_qc['OneReadPair']/pbc_qc['DistinctReadPairs']
pbc_qc['PBC2'] = pbc_qc['OneReadPair']/pbc_qc['TwoReadPairs']

pbc_qc.to_csv(path=pbc_qc_txt,index=True,sep='\\t')