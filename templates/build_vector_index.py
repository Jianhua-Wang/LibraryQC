#!/usr/bin/env python3

import sys
from subprocess import call

if len(sys.argv) == 1:
    sys.argv=["build_vector_index.py","$vector"]

txt = sys.argv[1]

fasta = 'vector.fasta'

with open(txt, 'r') as fp:
    ith = 0
    for line in fp:
        ith += 1
        if ith == 2:
            left = line[:-1]
        if ith == 4:
            right = line[:-1]

with open(fasta, 'w') as f:
    f.write('>vector\\n'+left+right)

call('bwa index vector.fasta', shell=True)