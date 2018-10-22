#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import re

if len(sys.argv) == 1:
    sys.argv=["extract_sg_from_SEfq.py","$fq1","$vector_txt","$left_shot","$right_shot"]


fastq_path = sys.argv[1]
vector_path = sys.argv[2]
left_shot = sys.argv[3]
right_shot = sys.argv[4]

txt_path = 'extract.txt'
failed_match_fq_path = 'fail_pattern_mapping.fastq'
matched_fq_path = 'succeed_pattern_mapping.fastq'


matched_raw = ['readID','sequence','distance','rc','quality']

def get_pattern(vector_path):
    with open(vector_path, 'r') as fp:
        ith = 0
        for line in fp:
            ith += 1
            if ith == 2:
                left = line[:-1]
            if ith == 4:
                right = line[:-1]

    left_pattern = left[-int(left_shot):].upper()
    right_pattern = right[:int(right_shot)].upper()

    pattern = '([ATCG]*%s)([ACTG]{20})(%s[ATCG]*)'%(left_pattern,right_pattern)
    return pattern

def DNA_rc(sequence):
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.replace('\\n', '')
    sequence = sequence.upper()
    return sequence[::-1]

def write_line(matched_raw,record,match):
    matched_raw[0] = record[0]
    matched_raw[1] = match.group(2)
    matched_raw[2] = len(match.group(1))
    matched_raw[3] = 'no_rc'
    matched_raw[4] = record[2][matched_raw[2]:matched_raw[2]+len(match.group(2))]
    matched_raw = list(map (str, matched_raw))
    return matched_raw

def pattern_mapping(record,pattern,txt,failed_match_fq):
    match = re.match(r'{}'.format(pattern), record[1])
    if match:
        txt.write('\\t'.join(write_line(matched_raw,record,match))+'\\n')
    else:
        failed_match_fq.write(record[0]+'\\n'+record[1]+'\\n+\\n'+record[2]+'\\n')

def extract_sgRNA(fastq_path,txt_path):
    fq = open(fastq_path,'r')
    txt =open(txt_path,'w')
    failed_match_fq = open(failed_match_fq_path,'w')
    ith = 0
    record = ['readID','sequence','quality']
    txt_line = ['readID','sequence','distance','rc','quality']
    txt.write('\\t'.join(txt_line)+'\\n')
    for line in fq:
        ith += 1
        if ith % 4 == 1:
            record[0] = line[:-1]
        elif ith % 4 == 2:
            record[1] = line[:-1]
        elif ith % 4 == 3:
            continue
        else:
            record[2] = line[:-1]
            pattern_mapping(record,pattern,txt,failed_match_fq)

pattern = get_pattern(vector_path)
extract_sgRNA(fastq_path,txt_path)

extract_df = pd.read_csv(txt_path,sep='\\t')
exact_distance = extract_df['distance'].mode()[0]

matched_fq = open(matched_fq_path, 'w')
for i in extract_df.index:
    line = extract_df.iloc[i].values
    matched_fq.write(line[0]+'\\n'+line[1]+'\\n+\\n'+line[4]+'\\n')