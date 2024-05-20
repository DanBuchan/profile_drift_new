import torch
from transformers import T5Tokenizer
from transformers.models.t5.modeling_t5 import T5ForConditionalGeneration
import re
import sys
import requests
from tqdm.auto import tqdm
import numpy as np
import csv
from collections import defaultdict
import pprint
import os

# tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
# model = T5ForConditionalGeneration.from_pretrained('Rostlab/prot_t5_xl_uniref50')
# # device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
# device = torch.device('cpu')

# model = model.to(device)
# model = model.eval()

# sequence_labels = "AFQFHEEHGEVCPANWQPGAKTIVANPQDSH"
# sequence_labels  = " ".join(sequence_labels)
# label_seq = tokenizer([sequence_labels], return_tensors="pt").input_ids

# seqs = ["AFQFHEEUUUVCPANWQPGUUUIVANPQDSH",
#         "UFQFHEEUGEVCPANWUPGAKTUVANPQDSU",
#         "AFUFHUEHGEVCPANWUPGAKUIVANPQDSH",]

# # print(model(input_ids=input_seq, labels=label_seq).loss)

# for i in range(0, 3):
#     seq  = " ".join(seqs[i])
#     seq = re.sub(r"[UZOB]", "<extra_id_0>", seq) 
#     input_seq = tokenizer([seq], return_tensors="pt").input_ids
#     output = model(input_ids=input_seq, labels=label_seq)
#     result = output["logits"].detach().numpy()
#     pred_array = np.argmax(result, axis=2)[0]
#     last_index = len(pred_array) -1
#     pred_array = np.delete(pred_array, last_index)
#     print(tokenizer.decode(pred_array))


# prottrans_seq_gen.py results_data/drift/drift_summary/alpha_fold_targets.csv

def get_pfam_seqs(pfam_fa, targets):
    print(targets)
    seqs = defaultdict(list)
    count = 0
    with open(pfam_fa, "r", encoding="utf-8") as fhIn:
        family_seqs = [] 
        current_family = ''
        header = ''
        seq = ''
        prt_ctl = False
        for line in fhIn:
            line = line.rstrip()
            if line.startswith(">"):
                this_family = line.split("|")[1]
                if prt_ctl:
                    # print(current_family, targets[0])
                    if current_family in targets:
                        seqs[current_family].append({
                            "header": header,
                            "seq": seq
                        })
                        count+=1
                        if count == 10:
                            pprint.pp(seqs)
                            exit()
                    header = ''
                
                current_family = this_family
                seq = ''
                header = line
                prt_ctl = True
            else:
                seq += line
        if current_family in targets:
            seqs[current_family].append({
                "header": header,
                "seq": seq
            })

    return seqs    

def read_targets(targets):
    data = set()
    with open(targets, "r", encoding="utf-8") as fhIn:
        targetreader = csv.reader(fhIn, delimiter=',')
        next(targetreader)
        for row in targetreader:
            data.add(row[0])
    return(data)

target_families = list(read_targets(sys.argv[1]))
# print(target_families)
if os.path.isfile('pfam_targets_for_prottrans.fa'):
    pfam_seqs = get_pfam_seqs('pfam_targets_for_prottrans.fa', target_families)
else:
    pfam_seqs = get_pfam_seqs(sys.argv[2], target_families)
    
    with open('pfam_targets_for_prottrans.fa', "w", encoding="utf-8") as fhOut:
        for family in pfam_seqs:
            print(len(pfam_seqs[family]))
            for seq_data in pfam_seqs[family]:
                header = seq_data['header']
                seq = seq_data['seq']
                fhOut.write(f'{header}\n{seq}\n')
