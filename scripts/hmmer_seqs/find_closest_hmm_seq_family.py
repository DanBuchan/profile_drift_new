import sys
from collections import defaultdict
import re
import pprint
import glob
import csv
from os.path import exists

# usage: python scripts/hmmer_seqs/find_closest_hmm_seq_family.py results_data/hmmer_matches/hmm_generated_seqs.fa results_data/psiblast_iteration_summaries/  ~/Data/pfam/Pfam-A.full.uniprot.fa

def get_hmm_generated_sequences(hmm_seqs):
    seqs = defaultdict(list)
    with open(hmm_seqs, "r", encoding="utf-8") as fhIn:
        header = ''
        seq = ''
        prt_ctl = False
        for line in fhIn:
            line = line.strip()
            # print(line)
            if line.startswith(">"):
                match = re.search("^>.+\|(PF\d+)-sample\d+", header)
                pf_family=''
                if match:
                    # print(match)
                    pf_family = match.groups()[0]
                if prt_ctl:
                    seqs[pf_family].append({"header": header,
                                            "seq": seq})
                seq = ''
                header = line
                prt_ctl = True
            else:
                seq += line 
    return seqs       


def get_pfam_family_names(seqs):
    family_names = set()
    for pf_id in seqs:
        for seq_data in seqs[pf_id]:
            match = re.search("^>(.+)\|PF\d+-sample\d+", seq_data['header'])
            if match:
                family_names.add(match.groups()[0])
    return list(family_names)


def collect_all_fasta_seqs(pfam_family_names, summaries, pfam_fasta):
    all_families = set()
    for family_name in pfam_family_names:
        for file in glob.glob(f'{summaries}/*{family_name}-*.csv'):
            with open(file, "r", encoding="utf-8") as fhIn:
                summary = csv.reader(fhIn, delimiter=',')
                for row in summary:
                    all_families.add(row[1])
                    all_families.add(row[2])
    # print(all_families)
    with open(pfam_fasta, "r", encoding="utf-8") as fhIn:
        header = ''
        seq = ''
        prt_ctl = False
        fhOut = open("all_drift_family_seqs.fa", "w", encoding="utf=8")
        for line in fhIn:
            line = line.strip()
            # print(line)
            if line.startswith(">"):
                match = re.search("^>.+\|(PF\d+)", header)
                pf_family=''
                if match:
                    pf_family = match.groups()[0]
                if prt_ctl and pf_family in all_families:
                    fhOut.write(f'{header}\n')
                    fhOut.write(f'{seq}\n')
                header = line
                seq = ''
                prt_ctl = True
            else:
                seq += line 


def read_fasta_db_seqs(fasta_file):
    all_seqs = defaultdict(list)
    with open(fasta_file, "r", encoding="utf-8") as fhIn:
        header = ''
        seq = ''
        prt_ctl = False
        for line in fhIn:
            line = line.strip()
            # print(line)
            if line.startswith(">"):
                match = re.search("^>.+\|(PF\d+)\d+", header)
                pf_family=''
                if match:
                    # print(match)
                    pf_family = match.groups()[0]
                if prt_ctl:
                    all_seqs[pf_family].append({"header": header,
                                            "seq": seq})
                seq = ''
                header = line
                prt_ctl = True
            else:
                seq += line
    return(all_seqs)


def find_closest_fasta(all_family_seqs, summaries, generated_seqs):
    for family in generated_seqs:
        first_entry = generated_seqs[family][0]
        print(first_entry)
        exit()

# 1. open file of generated seqs, read in and get family ID etc
generated_seqs = get_hmm_generated_sequences(sys.argv[1])
pprint.pp(generated_seqs)
# {PFID: [{header:,
#          seq: },]}
# 2. get a unique list of the pam family names
pfam_family_names = get_pfam_family_names(generated_seqs)
pprint.pp(generated_seqs)
# 3. Collect a fasta file of all the seqs we are going to need
if not exists("all_drift_family_seqs.fa"):
    collect_all_fasta_seqs(pfam_family_names, sys.argv[2], sys.argv[3])
all_family_seqs = read_fasta_db_seqs("all_drift_family_seqs.fa")
# 4. loop over the generated seqs, for each family make a fasta db. Then fasta all the seqs
find_closest_fasta(all_family_seqs, sys.argv[2], generated_seqs)

# pprint.pp(pfam_family_names)