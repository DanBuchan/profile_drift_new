import sys
from collections import defaultdict
import re
import pprint
import glob
import csv
from os.path import exists
from subprocess import Popen, PIPE

# usage: python scripts/hmmer_seqs/find_closest_hmm_seq_family.py results_data/hmmer_matches/hmm_generated_seqs.fa ~/Data/pfam/Pfam-A.full.uniprot.fa

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

def find_closest_fasta(all_family_seqs, generated_seqs):
    for pf_id in generated_seqs:        
        for seq_record in generated_seqs[pf_id]:
            match = re.search("^>.+\|(PF\d+-sample)\d+", seq_record['header'])
            query_name = ''
            if match:
                query_name = match.groups()[0]
            fhtmp = open(f"{query_name}.fa", "w", encoding="utf-8")
            fhtmp.write(f'{seq_record["header"]}\n')
            fhtmp.write(f'{seq_record["seq"]}\n')
            fhtmp.close()
            args = ['/home/dbuchan/Applications/fasta36/bin/fasta36',
                    '-q',
                    '-p',
                    '-O',
                    'out',
                    f"{query_name}.fa", 
                    f"{all_family_seqs}",
            ]
            print("Calculating", " ".join(args))
            try:
                p = Popen(args, stdout=PIPE, stderr=PIPE)
                result_stdout, err = p.communicate()
            except Exception as e:
                print(str(e))
                sys.exit(1)
            if p.returncode != 0:
                print("Non Zero Exit status: "+str(p.returncode))
                raise OSError("Non Zero Exit status: "+str(p.returncode))
            results = result_stdout.decode('utf-8')
            lines = results.split("\n")
            parse_results = False
            best_hit = ''
            best_score = ''
            for line in lines:
                if parse_results:
                    entries = line.split()
                    best_hit = entries[0]
                    try:
                        best_hit = best_hit.split("___")[1]
                        best_score = line[62:]
                        best_score = float(best_score.split()[0])
                    except:
                        best_hit = 'NA'
                        best_score = 'NA'
                    break
                if line.startswith("The best scores are:"):
                    parse_results = True
                if "residues in 1 query   sequences" in line:
                    parse_results = False
            print(query_name, best_hit, best_score)
            exit()

# 1. open file of generated seqs, read in and get family ID etc
generated_seqs = get_hmm_generated_sequences(sys.argv[1])
# pprint.pp(generated_seqs)
# {PFID: [{header:,
#          seq: },]}
# Header format ">pfam_family_name|PFID-sampleNN"
find_closest_fasta(sys.argv[2], generated_seqs)):