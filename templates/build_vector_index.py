#!/usr/bin/env python3

import sys
from subprocess import call

if len(sys.argv) == 1:
    sys.argv=["build_vector_index.py","$vector"]

txt = sys.argv[1]

fasta = 'vector.fasta'

with open(fasta, 'w') as f:
    with open(txt, 'r') as fp:
        for line in fp:
            f.write('>vector\\n' + line)

call('bowtie2-build -f vector.fasta vector', shell=True)