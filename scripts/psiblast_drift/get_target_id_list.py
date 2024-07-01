from collections import defaultdict
import sys
import csv

# python get_target_id_list.py ~/Data/pfam/pfam_consensus_reps_labelled_flattened.fa results_data/drift_summary/alpha_fold_targets.csv


def read_fasta(file):
    seqs = {}
    fasta_name = None
    seq = None
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                if fasta_name:
                    seqs[fasta_name[1:]] = seq
                fasta_name = line.rstrip()
                seq = ''
            else:
                seq = seq+line.rstrip()
    return(seqs)
 
def read_targets(file):
    targets = []
    with open(file) as csvfile:
        targetreader = csv.reader(csvfile, delimiter=',') 
        next(targetreader)
        for row in targetreader:
            targets.append(row[0])
    return(targets)

rep_file = sys.argv[1]
rep_seqs = read_fasta(rep_file)
id_list = list(rep_seqs.keys())
# print(id_list)
target_file = sys.argv[2]
target_list = read_targets(target_file)

for i, id in enumerate(id_list):
    pfam_id = id[-7:]
    if pfam_id in target_list:
       #print(id, pfam_id, i+1)
       print(i+1)
       