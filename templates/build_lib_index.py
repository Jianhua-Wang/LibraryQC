#!/usr/bin/env python3

import sys
from subprocess import call

if len(sys.argv) == 1:
    sys.argv=["build_lib_index.py","$library"]

file = sys.argv[1]

fasta = 'lib.fasta'

with open(fasta, 'w') as f:
    with open(file, 'r') as fp:
        ith = 0
        for line in fp:
            ith += 1
            line = line.split(',')
            if ith == 1:
                continue
            if line[1] == '+':
                f.write('>' + str(ith-1) + ':+\\n' + line[0] + '\\n')
            elif line[1] == '-':
                f.write('>' + str(ith-1) + ':-\\n' + line[0] + '\\n')
            else:
                f.write('>' + str(ith-1) + ':na\\n' + line[0]  + '\\n')

call('bwa index lib.fasta', shell=True)