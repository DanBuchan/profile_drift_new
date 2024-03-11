import sys
from subprocess import Popen, PIPE
import os
import csv

# USAGE
# python pfam_reps_blast.py ~/Data/pfam/pfam_consensus_reps.fa
#

def execute_process(executable_args, stdout_location=None):
    try:
        # print(' '.join(executable_args))
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
        fstdout = open(stdout_location, "w")
        fstdout.write(output.decode("utf-8") )
        fstdout.close()


def run_nw(rep_seqs):
    for seqa in rep_seqs.keys():
        for seqb in rep_seqs.keys():
            if seqa == seqb:
                continue
            fhOut = open('tmpa.fa', "w")
            fhOut.write(f'{seqa}\n')
            fhOut.write(f'{rep_seqs[seqa]}')
            fhOut.close()
            fhOut = open('tmpb.fa', "w")
            fhOut.write(f'{seqb}\n')
            fhOut.write(f'{rep_seqs[seqb]}')
            fhOut.close()
            nw_args = [
               '/home/dbuchan/Applications/EMBOSS-6.4.0/emboss/needle',
               'tmpa.fa',
               'tmpb.fa',
               'out.needle',
               '-gapopen=10',
               '-gapextend=0'
            ]
            execute_process(nw_args)
            os.remove('tmpa.fa')
            os.remove('tmpb.fa')
            with open('out.needle', "r") as fhIn:
                for line in fhIn:
                    if line.startswith('# Score:'):
                        print(f'{seqa[1:]},{seqb[1:]},{line.rstrip()[9:]}')
            os.remove(f'out.needle')
        exit()


reps_file = sys.argv[1]
rep_seqs = {}
fasta_name = None
seq = None
with open(reps_file, "r") as fhIn:
    for line in fhIn:
        if line.startswith(">"):
            if fasta_name:
                rep_seqs[fasta_name] = seq
            fasta_name = line.rstrip()
            seq = ''
        else:
            seq = seq+line.rstrip()
run_nw(rep_seqs)
