'''
We calculate the NW alignment with MAFFT then the distance
calculation with RAxML.

python3 ./calculate_pfam_distances.py ~/Data/pfam/pfam_consensus_reps_labelled.fa 
'''

import subprocess
from subprocess import Popen, PIPE
import sys
import os
import numpy as np
from collections import defaultdict
from time import time
from multiprocessing import Pool, Array
import itertools


def select_string(seqlist):
    length = 1000000
    ave = sum(map(len, seqlist)) / len(seqlist)
    selected = ''
    for seq in seqlist:
        if abs(len(seq)-ave) < length:
            length = abs(len(seq)-ave)
            selected = seq
    return(seq)

def execute_process(executable_args, stdout_location=None):
    try:
        print(' '.join(executable_args))
        p = Popen(executable_args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        code = p.returncode
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if code != 0:
        print("Non Zero Exit status: "+str(code))
        print("Command:"+' '.join(executable_args))
        if code != 255:
            raise OSError("Non Zero Exit status: "+str(code)+output.decode('utf-8'))
    if stdout_location:
        fhalign = open(stdout_location, "w")
        fhalign.write(output.decode("utf-8") )
        fhalign.close()

def process_distances(rep_file):
    rep_seqs = []

    align_file = 'reps.afa'
    dist_file = 'reps.dist'
    mafft_args = ['/usr/local/bin/mafft',
                   rep_file]
    execute_process(mafft_args, align_file)
    # run raxml

    raxml_args = ['/home/dbuchan/Applications/standard-RAxML/raxmlHPC-PTHREADS-AVX',
                  '-T',
                  '4',
                  '-s',
                  align_file,
                  '-n',
                  dist_file,
                  '-m',
                  'PROTGAMMABLOSUM62',
                  '-N'
                  '2',
                  '-p',
                  '123',
                  '-f',
                  'x']
    execute_process(raxml_args)

    try:
        os.remove(align_file+".reduced")
    except:
        pass
    try:
        os.remove(f'RAxML_info.rep.dist')
    except:
        pass
    try:
        os.remove(f'RAxML_parsimonyTree.rep.dist.RUN.0')
    except:
        pass


input_file = sys.argv[1]

process_distances(input_file)
