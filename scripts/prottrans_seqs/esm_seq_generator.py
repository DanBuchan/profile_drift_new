from typing import List, Tuple, Optional, Dict,  Union, Callable
import sys
import string
from pathlib import Path
import csv

import numpy as np
import torch
from scipy.spatial.distance import squareform, pdist, cdist
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import SeqIO
import biotite.structure as bs
from biotite.structure.io.pdbx import PDBxFile, get_structure
from biotite.database import rcsb
from tqdm import tqdm
import pandas as pd
from collections import defaultdict

import esm

"""
python esm_seq_generator.py ../iteration_summary.csv ~/Data/pfam/Pfam-A.full.uniprot > out 2> err &
"""

torch.set_grad_enabled(False)

# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)

def read_sequence(filename: str) -> Tuple[str, str]:
    """ Reads the first (reference) sequences from a fasta or MSA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return record.description, str(record.seq)

def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)

def read_msa(filename: str) -> List[Tuple[str, str]]:
    """ Reads the sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq))) for record in SeqIO.parse(filename, "fasta")]

def extend(a, b, c, L, A, D):
    """
    input:  3 coords (a,b,c), (L)ength, (A)ngle, and (D)ihedral
    output: 4th coord
    """

    def normalize(x):
        return x / np.linalg.norm(x, ord=2, axis=-1, keepdims=True)

    bc = normalize(b - c)
    n = normalize(np.cross(b - a, bc))
    m = [bc, np.cross(n, bc), n]
    d = [L * np.cos(A), L * np.sin(A) * np.cos(D), -L * np.sin(A) * np.sin(D)]
    return c + sum([m * d for m, d in zip(m, d)])

# Select sequences from the MSA to maximize the hamming distance
# Alternatively, can use hhfilter 
def greedy_select(msa: List[Tuple[str, str]], num_seqs: int, mode: str = "max") -> List[Tuple[str, str]]:
    assert mode in ("max", "min")
    if len(msa) <= num_seqs:
        return msa
    
    array = np.array([list(seq) for _, seq in msa], dtype=np.bytes_).view(np.uint8)

    optfunc = np.argmax if mode == "max" else np.argmin
    all_indices = np.arange(len(msa))
    indices = [0]
    pairwise_distances = np.zeros((0, len(msa)))
    for _ in range(num_seqs - 1):
        dist = cdist(array[indices[-1:]], array, "hamming")
        pairwise_distances = np.concatenate([pairwise_distances, dist])
        shifted_distance = np.delete(pairwise_distances, indices, axis=1).mean(0)
        shifted_index = optfunc(shifted_distance)
        index = np.delete(all_indices, indices)[shifted_index]
        indices.append(index)
    indices = sorted(indices)
    return [msa[idx] for idx in indices]

def compute_precisions(
    predictions: torch.Tensor,
    targets: torch.Tensor,
    src_lengths: Optional[torch.Tensor] = None,
    minsep: int = 6,
    maxsep: Optional[int] = None,
    override_length: Optional[int] = None,  # for casp
):
    if isinstance(predictions, np.ndarray):
        predictions = torch.from_numpy(predictions)
    if isinstance(targets, np.ndarray):
        targets = torch.from_numpy(targets)
    if predictions.dim() == 2:
        predictions = predictions.unsqueeze(0)
    if targets.dim() == 2:
        targets = targets.unsqueeze(0)
    override_length = (targets[0, 0] >= 0).sum()

    # Check sizes
    if predictions.size() != targets.size():
        raise ValueError(
            f"Size mismatch. Received predictions of size {predictions.size()}, "
            f"targets of size {targets.size()}"
        )
    device = predictions.device

    batch_size, seqlen, _ = predictions.size()
    seqlen_range = torch.arange(seqlen, device=device)

    sep = seqlen_range.unsqueeze(0) - seqlen_range.unsqueeze(1)
    sep = sep.unsqueeze(0)
    valid_mask = sep >= minsep
    valid_mask = valid_mask & (targets >= 0)  # negative targets are invalid

    if maxsep is not None:
        valid_mask &= sep < maxsep

    if src_lengths is not None:
        valid = seqlen_range.unsqueeze(0) < src_lengths.unsqueeze(1)
        valid_mask &= valid.unsqueeze(1) & valid.unsqueeze(2)
    else:
        src_lengths = torch.full([batch_size], seqlen, device=device, dtype=torch.long)

    predictions = predictions.masked_fill(~valid_mask, float("-inf"))

    x_ind, y_ind = np.triu_indices(seqlen, minsep)
    predictions_upper = predictions[:, x_ind, y_ind]
    targets_upper = targets[:, x_ind, y_ind]

    topk = seqlen if override_length is None else max(seqlen, override_length)
    indices = predictions_upper.argsort(dim=-1, descending=True)[:, :topk]
    topk_targets = targets_upper[torch.arange(batch_size).unsqueeze(1), indices]
    if topk_targets.size(1) < topk:
        topk_targets = F.pad(topk_targets, [0, topk - topk_targets.size(1)])

    cumulative_dist = topk_targets.type_as(predictions).cumsum(-1)

    gather_lengths = src_lengths.unsqueeze(1)
    if override_length is not None:
        gather_lengths = override_length * torch.ones_like(
            gather_lengths, device=device
        )

    gather_indices = (
        torch.arange(0.1, 1.1, 0.1, device=device).unsqueeze(0) * gather_lengths
    ).type(torch.long) - 1

    binned_cumulative_dist = cumulative_dist.gather(1, gather_indices)
    binned_precisions = binned_cumulative_dist / (gather_indices + 1).type_as(
        binned_cumulative_dist
    )

    pl5 = binned_precisions[:, 1]
    pl2 = binned_precisions[:, 4]
    pl = binned_precisions[:, 9]
    auc = binned_precisions.mean(-1)

    return {"AUC": auc, "P@L": pl, "P@L2": pl2, "P@L5": pl5}


def evaluate_prediction(
    predictions: torch.Tensor,
    targets: torch.Tensor,
) -> Dict[str, float]:
    if isinstance(targets, np.ndarray):
        targets = torch.from_numpy(targets)
    contact_ranges = [
        ("local", 3, 6),
        ("short", 6, 12),
        ("medium", 12, 24),
        ("long", 24, None),
    ]
    metrics = {}
    targets = targets.to(predictions.device)
    for name, minsep, maxsep in contact_ranges:
        rangemetrics = compute_precisions(
            predictions,
            targets,
            minsep=minsep,
            maxsep=maxsep,
        )
        for key, val in rangemetrics.items():
            metrics[f"{name}_{key}"] = val.item()
    return metrics

# read in Pfam stockholm data ~/data/pfam/Pfam-A.full.uniprot
def generate_seqs(msa, msa_transformer, msa_transformer_alphabet, align_name, mask_amount, fhOut):

    msa_transformer_batch_converter = msa_transformer_alphabet.get_batch_converter()
    msa_transformer_predictions = {}
    msa_transformer_results = []
    results = defaultdict(float)
    for name, inputs in msa.items():
        inputs = greedy_select(inputs, num_seqs=200) # can change this to pass more/fewer sequence
        msa_transformer_batch_labels, msa_transformer_batch_strs, msa_transformer_batch_tokens = msa_transformer_batch_converter([inputs])
        input_tokens = msa_transformer_batch_tokens.cpu().numpy()[0]
        if input_tokens.shape[1] > 1024:
            # we don't generate seqs longer than 1024 residues
            print(f"Skipping Familiy: {align_name}")
            continue
        substitution_numbers = round(len(input_tokens[0])*mask_amount)
        mask = torch.rand(msa_transformer_batch_tokens.shape).argsort(2) < substitution_numbers
        msa_transformer_batch_tokens = torch.where(mask, 31, msa_transformer_batch_tokens)
        # print("init", msa_transformer_batch_tokens)
        for test_seq in msa_transformer_batch_tokens[0]:
            # print(test_seq)
            test_seq[0] = 0
        # print("alte", msa_transformer_batch_tokens)
        # print(msa_transformer_batch_labels)
        # print(msa_transformer_batch_strs)
        # print(input_tokens)
        #print(msa_transformer_batch_tokens)
        # HERE MASK n(75%) tokens - 0 is the masking/padding token?
        # Char index (see data.py from_architerture()):
        #    <cls>    0
        #    <pad>    1
        #    <eos>    2
        #    <unk>    3
        #    residues in constants.py 4-30, '.' is 29, '-' is 30
        #    <mask>   31
        # Mask IDX is 31
        msa_transformer_batch_tokens = msa_transformer_batch_tokens.to(next(msa_transformer.parameters()).device)
        msa_transformer_predictions[name] = msa_transformer(msa_transformer_batch_tokens)
        # print(msa_transformer_predictions[name]['logits'])
        # print(msa_transformer_predictions[name]['logits'].size())
        input_tokens = msa_transformer_batch_tokens.cpu().numpy()[0]
        preds = msa_transformer_predictions[name]['logits'].cpu().numpy()
        for result in msa_transformer_predictions[name]['logits'].cpu().numpy():
            for i, seq in enumerate(result):
                # print("comparing")
                # print(input_tokens[i])
                pred_array = np.argmax(seq, axis=1)
                # print(pred_array)
                output_seq = ''
                for token in pred_array:
                    output_seq += msa_transformer_alphabet.get_tok(token)
                    output_seq = output_seq.replace(".", "")
                    output_seq = output_seq.replace("-", "")
                    output_seq = output_seq.replace("-", "")
                fhOut.write(f">{align_name}_{i}\n")
                fhOut.write(f"{output_seq}\n")
                # print(pred_array)
                tp_count = np.sum(input_tokens[i] == pred_array)
                pred_size = len(input_tokens[i])
                tpr = tp_count/pred_size
                # print(f'{name} {i} tpr: {tpr}: {pred_size}')
                results[name] += tpr
            results[name] = results[name]/(i+1)
            if not i == 199:
                print(f"Less than 200 seqs generated for: {align_name}, n = {i}")
    print(results)
   
def read_pfam_alignments(file, drift_families, msa_transformer, msa_transformer_alphabet):
    align_count = 0

    fh25 = open("masked_25_msa_transformer.fa", "w")
    fh50 = open("masked_50_msa_transformer.fa", "w")
    fh75 = open("masked_75_msa_transformer.fa", "w")
    with open(file, "rb") as fh:
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
                    if align_name in drift_families:
                        print(f"Processing: {align_name}")
                        # print(msa)
                        # print(len(msa))
                        for mask_amount in [[0.25, fh25], [0.5, fh50], [0.75, fh75]]:
                            generate_seqs(msa, msa_transformer, msa_transformer_alphabet, align_name, mask_amount[0], mask_amount[1])
                        fh25.flush()
                        fh50.flush()
                        fh75.flush()

                        # exit()
                    # run generator 
                    # reinitialise
                else:
                    align_count+=1
                align_name = ''
                msa = defaultdict(list)
            if line.startswith("#=GF AC   "):
                align_name = line[10:].rstrip()
                align_name = align_name.split(".", 1)[0]
            if not line.startswith("#"):
                entries = line.split()
                seq_data = (entries[0], remove_insertions(entries[1]))
                msa[align_name].append(seq_data)
    fh25.close()
    fh50.close()
    fh75.close()

def get_drift_set(file):
    drifts = set()
    with open(file, "r") as fhIn:
        next(fhIn)
        iteration_reader = csv.reader(fhIn, delimiter=",")
        for row in iteration_reader:
            drifts.add(row[3])
    return list(drifts)    


drift_families = get_drift_set(sys.argv[1])
msa_transformer, msa_transformer_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()

# remove contact head
# https://stackoverflow.com/questions/52548174/how-to-remove-the-last-fc-layer-from-a-resnet-model-in-pytorch
# layer_list = list(msa_transformer.children())
# layer_list.pop(3)
# msa_transformer = torch.nn.Sequential(*layer_list)
for name, param in msa_transformer.named_parameters():
    # print(name)
    param.requires_grad = False
msa_transformer = msa_transformer.eval().cuda()

msa_transformer_batch_converter = msa_transformer_alphabet.get_batch_converter()
# READ LIST OF FAMILIES THAT HAVE DRIFT

# python esm_seq_generator.py ../iteration_summary.csv ~/Data/pfam/Pfam-A.full.uniprot
read_pfam_alignments(sys.argv[2], drift_families, msa_transformer, msa_transformer_alphabet)