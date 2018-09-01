#!/usr/bin/env python3

import sys
import re

if len(sys.argv) == 1:
    sys.argv=["extract_sg_from_sam.py","$sam","$pattern_str"]

sam = sys.argv[1]
pattern = sys.argv[2]

txt = 'extract.txt'
failed_match_fq = 'fail_pattern_mapping.fastq'
matched_fq = 'succeed_pattern_mapping.fastq'

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
    matched_raw[4] = DNA_rc(read[10][::-1])[matched_raw[2]:matched_raw[2]+20]
    matched_raw = list(map (str, matched_raw))
    return matched_raw

def pattern_mapping_2reads(read1,read2):
    read1 = read1.split()
    read2 = read2.split()
    match1 = re.match(r'{}'.format(pattern), DNA_rc(read1[9]))
    match2 = re.match(r'{}'.format(pattern), DNA_rc(read2[9]))
    txt_line = ''
    fq_line = ''
    dup_line = ''
    if match1 and not match2:
        txt_line = '\\t'.join(write_line(matched_raw,read1,match1))+'\\n'
    elif match2 and not match1:
        txt_line = '\\t'.join(write_line(matched_raw,read2,match2))+'\\n'
    elif not match1 and not match2:
        fq_line = '@'+read1[0]+' 1:N:0:ACAAGCTA\\n'+read1[9]+'\\n+\\n'+read1[10]+'\\n@'+read2[0]+' 2:N:0:ACAAGCTA\\n'+read2[9]+'\\n+\\n'+read2[10]+'\\n'
    else:
        txt_line = '\\t'.join(write_line(matched_raw,read1,match1))+'\\n'+'\\t'.join(write_line(matched_raw,read2,match2))+'\\n'
    return txt_line,fq_line

failed_match = open(failed_match_fq, 'w')
extract_txt = open(txt, 'w')
with open(sam, 'r') as fp:
    extract_txt.write('\\t'.join(matched_raw)+'\\n')
    ith = 0
    matched_num = 0
    failed_num = 0
    for line in fp:
        if  re.match(r'^@',line):
            continue
        ith += 1
        if ith % 2 == 1:
            read1 = line
        if ith % 2 == 0:
            read2 = line
            mapping = pattern_mapping_2reads(read1,read2)
            if mapping[0] != '':
                extract_txt.write(mapping[0])
            else:
                failed_match.write(mapping[1])

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