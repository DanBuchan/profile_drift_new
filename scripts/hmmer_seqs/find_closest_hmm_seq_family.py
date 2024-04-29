import sys
from collections import defaultdict
import re
import pprint
import glob

# 2. Open MSA summary for family and get all families
# 3. make fasta db of all the families
# 4. Iteratr over seqs in family seeing what the top hits are

def get_hmm_generated_sequnces(hmm_seqs):
    seqs = defaultdict(list)
    with open(hmm_seqs, "r", encoding="utf-8") as fhIn:
        header = ''
        seq = ''
        prt_ctl = False
        for line in fhIn:
            line = line.strip()
            # print(line)
            if line.startswith(">"):
                match = re.search("^>.+\|(PF\d+)-sample\d+", line)
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

# 1. open file of generated seqs, read in and get family ID etc
generated_seqs = get_hmm_generated_sequnces(sys.argv[1])
# 2. get a unique list of the pfam family names
pfam_family_names = get_pfam_family_names(generated_seqs)
pprint.pp(pfam_family_names)