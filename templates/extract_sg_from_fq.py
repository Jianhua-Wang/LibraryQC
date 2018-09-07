#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import re

if len(sys.argv) == 1:
    sys.argv=["extract_sg_from_sam.py","$fq1","$fq2","$pattern_str"]

fq1 = sys.argv[1]
fq2 = sys.argv[2]
pattern = sys.argv[3]

two_in_one = 'two_in_one.txt'
txt = 'extract.txt'
failed_match_fq = 'fail_pattern_mapping.fastq'
matched_fq = 'succeed_pattern_mapping.fastq'

two_in_one_txt = open(two_in_one, 'w')

with open(fq1,'r') as f1, open(fq2,'r') as f2:
    ith = 0
    for x, y in zip(f1, f2):
        ith += 1
        if ith == 1:
            two_in_one_txt.write('readid'+'\\t'+'read1'+'\\t'+'quality1'+'\\t'+'read2'+'\\t'+'quality2'+'\\n')
            read_tail = x.split()[1][1:]
        if ith % 4 == 1:
            readid = x.split()[0][1:]
        if ith % 4 == 2:
            read1 = x[:-1]
            read2 = y[:-1]
        if ith % 4 == 0:
            quality1 = x[:-1]
            quality2 = y[:-1]
            two_in_one_txt.write(readid+'\\t'+read1+'\\t'+quality1+'\\t'+read2+'\\t'+quality2+'\\n')

two_in_one_txt.close()

def DNA_rc(sequence):
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.replace('\\n', '')
    sequence = sequence.upper()
    return sequence[::-1]

matched_raw = ['readID','sequence','distance','rc','quality']
failed_match_raw = ['readID','sequence','quality']

def write_line(matched_raw,read,match):
    matched_raw[0] = read[0]
    matched_raw[1] = match.group(2)
    matched_raw[2] = len(match.group(1))
    matched_raw[3] = 'rc'
    matched_raw[4] = DNA_rc(read[2])[matched_raw[2]:matched_raw[2]+20]
    matched_raw = list(map (str, matched_raw))
    return matched_raw

def pattern_mapping_2reads(read1,read2):
    match1 = re.match(r'{}'.format(pattern), DNA_rc(read1[1]))
    match2 = re.match(r'{}'.format(pattern), DNA_rc(read2[1]))
    txt_line = ''
    fq_line = ''
    if match1 and not match2:
        txt_line = '\\t'.join(write_line(matched_raw,read1,match1))+'\\n'
    elif match2 and not match1:
        txt_line = '\\t'.join(write_line(matched_raw,read2,match2))+'\\n'
    elif not match1 and not match2:
        fq_line = '@'+read1[0]+' 1'+read_tail+'\\n'+read1[1]+'\\n+\\n'+read1[2]+'\\n@'+read2[0]+' 2'+read_tail+'\\n'+read2[1]+'\\n+\\n'+read2[2]+'\\n'
    else:
        txt_line = '\\t'.join(write_line(matched_raw,read1,match1))+'\\n'+'\\t'.join(write_line(matched_raw,read2,match2))+'\\n'
    return txt_line,fq_line

failed_match = open(failed_match_fq, 'w')
extract_txt = open(txt, 'w')
with open(two_in_one, 'r') as fp:
    extract_txt.write('\\t'.join(matched_raw)+'\\n')
    ith = 0
    matched_num = 0
    failed_num = 0
    for line in fp:
        ith += 1
        if  ith == 1:
            continue
        line = line.split()
        read1 = [line[0],line[1],line[2]]
        read2 = [line[0],line[3],line[4]]
        mapping = pattern_mapping_2reads(read1,read2)
        if mapping[0] != '':
            extract_txt.write(mapping[0])
        else:
            failed_match.write(mapping[1])

df = pd.read_table(txt)
exact_distance = df['distance'].mode()[0]

fq = open(matched_fq, 'w')
with open(txt, 'r') as fp:
    ith = 0
    for line in fp:
        ith += 1
        if ith == 1:
            continue
        line = line.split()
        if line[2] == '58':
            fq.write('@'+line[0]+'\\n'+line[1]+'\\n+\\n'+line[4]+'\\n')

os.remove('two_in_one.txt')