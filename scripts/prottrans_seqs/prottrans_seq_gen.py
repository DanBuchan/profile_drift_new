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
import random
import math

# python prottrans_seq_gen.py ../../results_data/generation_or_af_targets/alphafold_targets.csv ~/Data/pfam/Pfam-A.full.uniprot.fa

def get_pfam_seqs(pfam_fa, targets):
    # print(targets)
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
                        # if count == 10:
                        #     pprint.pp(seqs)
                        #     exit()
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
            for seq_data in pfam_seqs[family]:
                header = seq_data['header']
                seq = seq_data['seq']
                fhOut.write(f'{header}\n{seq}\n')

def predict_seq(seq_labels, input_labels, t5, device):
    seq_labels  = " ".join(seq_labels)
    label_seq = tokenizer([seq_labels], return_tensors="pt").input_ids
    seq  = " ".join(input_labels)
    seq = re.sub(r"[UZOB]", "<extra_id_0>", seq) 
    input_seq = tokenizer([seq], return_tensors="pt").input_ids
    inputs = input_seq.to(device)
    labels = label_seq.to(device)
    output = t5(input_ids=inputs, labels=labels)
    del inputs
    del labels
    result = output["logits"].detach().cpu().numpy()
    pred_array = np.argmax(result, axis=2)[0]
    last_index = len(pred_array) -1
    pred_array = np.delete(pred_array, last_index)
    pred_str = "".join(list(tokenizer.decode(pred_array)))
    pred_str = pred_str.replace(" ", "")
    return pred_str

tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
model = T5ForConditionalGeneration.from_pretrained('Rostlab/prot_t5_xl_uniref50')
# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
device = torch.device('cuda')
model = model.to(device)
model = model.eval()
print("generating seqs")
masked_25 = open("masked_25_percent_targets.fa", "w", encoding="utf-8")
masked_50 = open("masked_50_percent_targets.fa", "w", encoding="utf-8")
masked_75 = open("masked_75_percent_targets.fa", "w", encoding="utf-8")
errors = open("dropped_targets.fa", "w", encoding="utf-8")

for family in pfam_seqs:
    print(family)
    # if("PF05086" not in family):
    #     continue
    sample = random.choices(pfam_seqs[family], k=25)
    # print(len(sample))
    for cnt, seq_data in enumerate(sample):
        torch.cuda.empty_cache()
        # print(seq_data)
        # print(f"OG     SEQ: {seq_data['seq']}")
        for i in [25, 50, 75]:
            fraction = i/100
            number = int(math.ceil(len(seq_data['seq']))*fraction)
            # print(number)
            locations = list(range(0, len(seq_data['seq'])))
            location_sample = random.sample(locations, k=number)
            new_seq = seq_data['seq']
            for location in location_sample:
                new_seq = new_seq[:location] + "U" + new_seq[location + 1:]
            
            # print(f"MASKED SEQ: {new_seq}")
            try:
                predicted_seq = predict_seq(seq_data['seq'], new_seq, model, device)
                predicted_seq = re.sub('<.+?>', '', predicted_seq)

                # print(f"PREDED SEQ: {predicted_seq}")
                if i == 25:
                    head = seq_data['header']
                    masked_25.write(f'{head}_masked{i}_{cnt}\n{predicted_seq}\n')
                if i == 50:
                    head = seq_data['header']
                    masked_50.write(f'{head}_masked{i}_{cnt}\n{predicted_seq}\n')
                if i == 75:
                    head = seq_data['header']
                    masked_75.write(f'{head}_masked{i}_{cnt}\n{predicted_seq}\n')
            except Exception as e:
                print(f"COULDN'T ALLOCATE {family}")
                errors.write(f"{seq_data['header']}\n")
                errors.write(f"{seq_data['seq']}\n")
            
        # exit()

masked_25.close()
masked_50.close()
masked_75.close()

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
