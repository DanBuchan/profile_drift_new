import sys
from subprocess import Popen, PIPE
import os
import csv

# USAGE
# python pfam_reps_nw.py ~/Data/pfam/pfam_consensus_reps.fa ~/Data/pfam/1_pfam_consensus
# Where [0] is the database to search
# And [1] is a file of fasta seqs to match to the db

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


def run_nw(dbseqs, repseqs):
    for seqa in repseqs.keys():
        for seqb in dbseqs.keys():
            if seqa == seqb:
                continue
            fhOut = open('tmpa.fa', "w")
            fhOut.write(f'{seqa}\n')
            fhOut.write(f'{dbseqs[seqa]}')
            fhOut.close()
            fhOut = open('tmpb.fa', "w")
            fhOut.write(f'{seqb}\n')
            fhOut.write(f'{repseqs[seqb]}')
            fhOut.close()
            nw_args = [
               '/home/dbuchan/EMBOSS-6.4.0/emboss/needle',
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


def read_fasta(file):
    seqs = {}
    fasta_name = None
    seq = None
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                if fasta_name:
                    seqs[fasta_name] = seq
                fasta_name = line.rstrip()
                seq = ''
            else:
                seq = seq+line.rstrip()
    return(seqs)

db_seqs = read_fasta(sys.argv[1])
rep_seqs = read_fasta(sys.argv[2])
run_nw(db_seqs, rep_seqs)
