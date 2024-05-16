import torch
from transformers import T5Tokenizer
from transformers.models.t5.modeling_t5 import T5ForConditionalGeneration
import re
import os
import requests
from tqdm.auto import tqdm
import numpy as np


tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
model = T5ForConditionalGeneration.from_pretrained('Rostlab/prot_t5_xl_uniref50')
# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
device = torch.device('cpu')

model = model.to(device)
model = model.eval()

sequence_labels = "AFQFHEEHGEVCPANWQPGAKTIVANPQDSH"
sequence_labels  = " ".join(sequence_labels)
label_seq = tokenizer([sequence_labels], return_tensors="pt").input_ids

seqs = ["AFQFHEEUUUVCPANWQPGUUUIVANPQDSH",
        "UFQFHEEUGEVCPANWUPGAKTUVANPQDSU",
        "AFUFHUEHGEVCPANWUPGAKUIVANPQDSH",]

# print(model(input_ids=input_seq, labels=label_seq).loss)

for i in range(0, 3):
    seq  = " ".join(seqs[i])
    seq = re.sub(r"[UZOB]", "<extra_id_0>", seq) 
    input_seq = tokenizer([seq], return_tensors="pt").input_ids
    output = model(input_ids=input_seq, labels=label_seq)
    result = output["logits"].detach().numpy()
    pred_array = np.argmax(result, axis=2)[0]
    last_index = len(pred_array) -1
    pred_array = np.delete(pred_array, last_index)
    print(tokenizer.decode(pred_array))


# https://github.com/agemagician/ProtTrans/issues/137