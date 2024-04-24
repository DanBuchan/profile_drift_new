from collections import defaultdict
import sys

# python run_raxml.py ~/Data/pfam/pfam_consensus_reps_labelled_flattened.fa results_data/drift_summary/alpha_fold_targets.csv


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

def read_targets(file):
    tagets = []
    return(targets)


rep_file = sys.argv[1]
rep_seqs = read_fasta(rep_file)
id_list = list(rep_seqs.keys())
print(id_list)
target_file = sys.argv[2]
read_targets(target_file)