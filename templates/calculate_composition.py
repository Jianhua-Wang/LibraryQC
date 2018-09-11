#!/usr/bin/env python3

import sys
import re
import numpy as np
import pandas as pd
import subprocess as sp

if len(sys.argv) == 1:
    sys.argv=["extract_fail_align_to_vector.py","$mismatch","$fastq","$extract","$right_position","$maptovector","$lib_csv"]

mismatch = sys.argv[1]
fastq = sys.argv[2]
extract = sys.argv[3]
right_position = sys.argv[4]
maptovector = sys.argv[5]
lib_csv = sys.argv[6]

ngs_composition_txt = 'ngs_composition.txt'
lib_composition_txt = 'lib_composition.txt'

def wc(file):
    return int(sp.check_output(["wc", "-l", file]).split()[0])

total = wc(fastq)/4
within_sgrna = wc(extract) - 1
exact = wc(right_position)/4
not_exact = within_sgrna - exact

mis_df = pd.read_table(mismatch)
mis_count = mis_df.groupby('mis').size()
mis0_lib_id = len(mis_df[mis_df.mis == '0'].groupby('lib_id').size())

with open(maptovector,'r') as fp:
    ith = 0
    empty = 0
    oneread_unmapped = 0
    tworead_unmapped = 0
    cigar = []
    readid = ''
    for line in fp:
        if  re.match(r'^@',line):
            continue
        ith += 1
        line_sp = line.split()
        if ith % 2 == 1:
            if line_sp[0] == readid:
                ith = ith - 1
                continue
            else:
                readid = line_sp[0]
                cigar.append(line_sp[5])
        else:
            cigar.append(line_sp[5])
            if len(cigar) > 2:
                print(ith)
            if '*' not in cigar:
                empty += 1
            elif cigar == ['*','*']:
                tworead_unmapped += 1
            else:
                oneread_unmapped += 1
            cigar = []

def calculate_per(a,total):
    return [a, str(a/total*100)[0:5]+'%']

ngs_composition = pd.DataFrame(columns=['Read count', 'Percentage'])

ngs_composition.loc['perfect map'] = calculate_per(mis_count[0],total)
ngs_composition.loc['mismatch1'] = calculate_per(mis_count[1],total)
ngs_composition.loc['mismatch2'] = calculate_per(mis_count[2],total)
ngs_composition.loc['unmapped'] = calculate_per(mis_count['unmapped'],total)
ngs_composition.loc['wrong position'] = calculate_per(not_exact,total)
ngs_composition.loc['empty'] = calculate_per(empty,total)
ngs_composition.loc['oneread_unmapped'] = calculate_per(oneread_unmapped,total)
ngs_composition.loc['tworead_unmapped'] = calculate_per(tworead_unmapped,total)
ngs_composition.loc['total'] = calculate_per(total,total)

ngs_composition.to_csv(ngs_composition_txt,index=True,sep='\\t')

total_lib = wc(lib_csv) - 1

def count_mis(a):
    return len(mis_df[mis_df.mis == a].groupby('lib_id').size())

lib_composition = pd.DataFrame(columns=['Lib count', 'Percentage'])
lib_composition.loc['perfect map'] = calculate_per(count_mis('0'),total_lib)
lib_composition.loc['mismatch1'] = calculate_per(count_mis('1'),total_lib)
lib_composition.loc['mismatch2'] = calculate_per(count_mis('2'),total_lib)
lib_composition.loc['unmapped'] = calculate_per(count_mis('unmapped'),total_lib)
lib_composition.loc['total_lib'] = calculate_per(total_lib,total_lib)

lib_composition.to_csv(lib_composition_txt,index=True,sep='\\t')