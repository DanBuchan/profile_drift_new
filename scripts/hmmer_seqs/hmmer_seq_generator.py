import csv
import sys
from collections import defaultdict

"""
python hmmer_seq_generator.py ../iteration_summary.csv ~/Data/pfam/Pfam-A.hmm
"""

def read_drifts(file):
    drift_families = defaultdict(set)
    with open(file, "r") as fhIn:
        next(fhIn)
        iteration_reader = csv.reader(fhIn, delimiter=",")
        for row in iteration_reader:
            drift_families[row[1]].add(row[3])
    return drift_families


def make_drift_set_non_redundant(drift_families):
    families = set()
    for key, hit_set in drift_families.items():
        families.add(key)
        families.update(hit_set)
    return(families)


def make_hmms_file(drift_set, hmm_file):
    current_hmm_data = []
    current_pfam_id = ''
    current_hmm_name = ''
    fhOut = open("hmm_subset.hmm", "w")
    with open(hmm_file, "r") as fhIn:
        for line in fhIn:
            current_hmm_data.append(line)
            if line.startswith("NAME"):
                current_hmm_name = line[6:]
            if line.startswith("ACC"):
                current_pfam_id = line[6:].split(".")[0]
            if line.startswith("//"):
                if current_pfam_id in drift_set:
                    replacement_name = current_hmm_data[1].rstrip()
                    replacement_name += f"|{current_pfam_id}\n"
                    current_hmm_data[1] = replacement_name
                    for output_ln in current_hmm_data:
                        fhOut.write(output_ln)
                current_hmm_data = []
                current_pfam_id = ''
    fhOut.close()

def runhmmemit():
    pass
    #NOT IMPLEMENTED, ran hmmemit on commandline as:
    # ~/Applications/hmmer-3.3.2/src/hmmemit -N 100 hmm_subset.hmm > hmm_generated_seqs.fa


drift_families = read_drifts(sys.argv[1])
nr_drift_set = make_drift_set_non_redundant(drift_families)
make_hmms_file(nr_drift_set, sys.argv[2])
# run_hmmemit()
