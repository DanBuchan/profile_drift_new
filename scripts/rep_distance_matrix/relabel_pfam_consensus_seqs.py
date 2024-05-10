import sys

'''
usage:
relabel_pfam_consensus_seqs.py ~/Data/pfam/Pfam-A.hmm ~/Data/pfam/pfam_consensus_reps.fa > pfam_consensus_reps_labelled.fa
'''

def read_pfam_hmm(input_file):
    names = dict()
    current_name = ''
    curren_acc = ''
    with open(input_file, "r") as fhIn:
        for line in fhIn:
            if line.startswith("NAME  "):
                current_name = line.rstrip()[6:]
            if line.startswith("ACC   "):
                current_acc= line.rstrip()[6:13]
                names[current_name] = current_acc
    return names

def annotate_fasta(names, fasta_file):

    print_ctl = False
    with open(fasta_file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                current_id = line[1:-11]
                if current_id in names:
                    print_ctl = True
                    print(f'{line[:-1]}|{names[current_id]}')
                else:
                    print_ctl = False
            else:
                if print_ctl:
                    print(line.rstrip())

pfam_names = read_pfam_hmm(sys.argv[1])
annotate_fasta(pfam_names, sys.argv[2])1-cysPrx_C-consensus|PF10417