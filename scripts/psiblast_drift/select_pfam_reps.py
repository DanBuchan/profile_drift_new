import sys
import random

# python select_pfam_reps.py ~/Data/pfam/Pfam-A.full.uniprot.fa

pfam_fa = sys.argv[1]

with open(pfam_fa, "r", encoding="utf-8") as fhIn:
    family_seqs = [] 
    current_family = ''
    header = ''
    seq = ''
    prt_ctl = False
    
    for line in fhIn:
        # print(line)
        line = line.rstrip()
        if line.startswith(">"):
            
            this_family = line.split("|")[1]
            if prt_ctl:
                family_seqs.append({'header': header, 'seq': seq})
                if this_family not in current_family:
                    current_family = this_family
                    ## print(family_seqs)
                    rand_rep = random.choice(family_seqs)
                    # print(rand_rep)
                    print(rand_rep['header'])
                    print(rand_rep['seq'])
                    family_seqs = []
                    ##exit()
            else:
                current_family = this_family
            seq = ''
            header = line
            prt_ctl = True
        else:
            seq += line

    family_seqs.append({'header': header, 'seq': seq})
    rand_rep = random.choice(family_seqs)
    print(rand_rep['header'])
    print(rand_rep['seq'])
    