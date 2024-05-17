import sys
import csv
from collections import defaultdict
from os.path import exists
import os
from subprocess import Popen, PIPE

"""
python find_closest_family_rep.py ../iteration_summary.csv ~/Data/pfam/Pfam-A.full.uniprot
"""

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_pfam_alignments(pfam_aligns, drift_families):
    # read through the PFAM A alignments and put them all in a directory
    nr_list = set()
    for family in drift_families:
         nr_list.add(family)
         for item in drift_families[family]:
             nr_list.add(item)

    align_count = 0
    with open(pfam_aligns, "rb") as fh:
        align_name = ''
        msa = defaultdict(list)
        for line_binary in fh:
            try:
                line = line_binary.decode("utf-8")
            except Exception as e:
                print(align_name, e)
                print(line_binary)
                continue
            if line.startswith("//"):
                continue
            if line.startswith("# STOCKHOLM"):
                if align_count != 0:
                    if align_name in nr_list:
                        print(f"Printing: {align_name}")
                        with open(f"{align_name}.fa", "w") as fhOut:
                            for msa_line in msa[align_name]:
                                fhOut.write(f">{msa_line[0]}\n")
                                fhOut.write(f"{msa_line[1]}\n")             
                else:
                    align_count+=1
                align_name = ''
                msa = defaultdict(list)
            if line.startswith("#=GF AC   "):
                align_name = line[10:].rstrip()
                align_name = align_name.split(".", 1)[0]
            if not line.startswith("#"):
                entries = line.split()
                seq = entries[1].replace('-', '')
                seq = seq.replace('.', '')
                seq_data = (entries[0], seq)
                msa[align_name].append(seq_data)
        print(f"Printing: {align_name}")
        with open(f"{align_name}.fa", "w") as fhOut:
            for msa_line in msa[align_name]:
                fhOut.write(f">{msa_line[0]}\n")
                fhOut.write(f"{msa_line[1]}\n")

    with open("families_list.txt", "w") as fhOut:
        for entry in nr_list:
            fhOut.write(f'{entry}\n')
        


def read_drifts(file):
    drift_families = defaultdict(set)
    with open(file, "r") as fhIn:
        next(fhIn)
        iteration_reader = csv.reader(fhIn, delimiter=",")
        for row in iteration_reader:
            drift_families[row[1]].add(row[3])
    return drift_families

def read_generated_seqs(file):
    seqs = defaultdict(list)
    family_id = ''
    current_prot_id = ''
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                current_prot_id = line[1:].strip()
                family_id = line[1:]
                family_id = family_id.rstrip()
                family_id = family_id.split("_")[0]
            else:
                seqs[family_id].append([current_prot_id, line.rstrip().replace("-", '')])
    return seqs

def read_generated_seqs_hmmer(file):
    seqs = defaultdict(list)
    family_id = ''
    current_prot_id = ''
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                current_prot_id = line[1:].strip()
                family_id = line[1:]
                family_id = family_id.rstrip()
                family_id = family_id.split("|")[1]
                family_id = family_id.split("-")[0]
            else:
                seqs[family_id].append([current_prot_id, line.rstrip().replace("-", '')])
    return seqs

def read_fasta_seqs(family_id, file):
    seqs = []
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                pass
            else:
                seqs.append(line.rstrip().replace("-", ''))
    return seqs
     
def run_fasta(seq_dyad):
    seq_name = seq_dyad[0]
    seq = seq_dyad[1]
    if len(seq) <= 1:
        return[seq_name, "NA", 0]
    with open("query.fa", "w") as fhOut: 
        fhOut.write(">{seq_name}\n")
        fhOut.write(f"{seq}\n")
    # for pair in pairs:
    args = ['/home/dbuchan/Applications/fasta36/bin/fasta36',
            '-q',
            '-p',
            '-O',
            'out',
            'query.fa', 
            f'./search_targets.fa'
            ]
    print("Calculating", " ".join(args))
    try:
        p = Popen(args, stdout=PIPE, stderr=PIPE)
        results, err = p.communicate()
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if p.returncode != 0:
        print("Non Zero Exit status: "+str(p.returncode))
        raise OSError("Non Zero Exit status: "+str(p.returncode))
    results = results.decode('utf-8')
    parse_results = False
    lines = results.split("\n")
    best_hit = ''
    best_score = ''
    for line in lines:
        if parse_results:
            entries = line.split()
            best_hit = entries[0]
            best_hit = best_hit.split("___")[1]
            best_score = line[62:]
            best_score = float(best_score.split()[0])
            break
        if line.startswith("The best scores are:"):
            parse_results = True
        if "residues in 1 query   sequences" in line:
            parse_results = False
    
    return [seq_name, best_hit, best_score]

# loop over every 
def find_closest_fasta(generated_seqs, alignment_list, pfam_family, families_hit):
    # target_seqs = {}
    proceed_analysis = True
    fhOut = open("search_targets.fa", "w")
    search_seqs = []
    for target in families_hit:
        # here we generate a searchable file.
        if target in alignment_list:
            seq_id = ''
            with open(f"alignments/{target}.fa") as fhIn:
                for line in fhIn:
                    if line.startswith(">"):
                        seq_id = line.rstrip()
                        seq_id = f"{seq_id}___{target}"
                    else:
                        search_seqs.append([seq_id,line.rstrip()])
            pass
            # target_seqs[target] = read_fasta_seqs(target, f"alignments/{target}.fa")
        else:
            proceed_analysis = False
            eprint(f"NOT ANALYSING {pfam_family} vs {target}")
    for seq in search_seqs:
        fhOut.write(f"{seq[0]}\n")
        fhOut.write(f"{seq[1]}\n")  
    fhOut.close()

    results = []
    if proceed_analysis:
        for seq in generated_seqs[pfam_family]:
            print("running fasta")
            best_hit = run_fasta(seq)
            results.append([pfam_family] + best_hit )
    return results

drift_families = read_drifts(sys.argv[1])
if not exists("families_list.txt"):
    parse_pfam_alignments(sys.argv[2], drift_families)

alignment_list = []
for file in os.listdir("./alignments"):
    if file.endswith(".fa"):
        alignment_list.append(file[:-3])


fhResults = open("summarised_msa_model_results.csv", "w")
fhResults.write("file,generated_family,query_name,best_hit_family,best_hit_score\n")
for file in ['hmm_generated_seqs_flattened.fa']:
    # generated_seqs = read_generated_seqs(file)
    generated_seqs = read_generated_seqs_hmmer(file)
    # print(generated_seqs)
    for pf_family in drift_families:
        # print(pf_family, drift_families[pf_family])
        results = find_closest_fasta(generated_seqs, alignment_list, pf_family, drift_families[pf_family])
        # print(results)
        for hit in results:
            # print(hit)
            fhResults.write(f"{file},{hit[0]},{hit[1]},{hit[2]},{hit[3]}\n")
            fhResults.flush()
        # exit()
