import sys
from collections import defaultdict
import re

def get_hmm_generated_sequnces(hmm_seqs):
    
    seqs = defaultdict(list)
    with open(hmm_seqs, "r", encoding="utf-8") as fhIn:
        header = ''
        seq = ''
        prt_ctl = False
        for line in fhIn:
            line = line.strip()
            print(line)
            if line.startswith(">"):
                match = re.search("^>.+\|(PF\d+)-sample\d+", line)
                if match:
                    print(match)
                if prt_ctl:
                    seqs[pf_family].append({"header": header,
                                            "seq": seq})
                    seq = ''
                header = line
                prt_ctl = True
            else:
                seq += line        

get_hmm_generated_sequnces(sys.argv[1])
