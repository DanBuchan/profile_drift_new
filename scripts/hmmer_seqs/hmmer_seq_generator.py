import sys
import glob
import csv
import subprocess

# python ./scripts/hmmer_seqs/hmmer_seq_generator.py /home/dbuchan/Projects/profile_drift/results_data/generation_or_af_targets/alphafold_targets_as_id.txt /home/dbuchan/Projects/profile_drift/results_data/alphafold_targets ~/Data/pfam/Pfam-A.hmm 50

# Read in list of targets alpha_fold_targets.csv
# Open drift summaries for each and get full list of PFIDs

def read_mafft_targets(file):
    values = []
    with open(file, "r") as fhIn:
        for line in fhIn:
            values.append(line.rstrip())
    return values

def read_summaries(targets, summaries):
    target_pfam_ids = set()
    all = []
    for target in targets:
        print(target)
        file = list(glob.glob(f"{summaries}/{target}_*.csv"))[0]
        # print(file)
        with open(file, "r", encoding="utf-8") as fhIn:
            next(fhIn)
            try:
                csvreader = csv.reader(fhIn, delimiter=',')
                entries = next(csvreader)
                target_pfam_ids.add(entries[1])
                all.append(entries[1])
                # print(entries)
            except Exception as e:
                pass
    # print(all)
    return(list(target_pfam_ids))


def make_hmms_file(drift_set, hmm_file):
    current_hmm_data = []
    current_pfam_id = ''
    current_hmm_name = ''
    fhOut = open("hmm_subset.hmm", "w")
    with open(hmm_file, "r") as fhIn:
        for line in fhIn:
            current_hmm_data.append(line)
            if line.startswith("NAME"):
                current_hmm_name = line[6:]
            if line.startswith("ACC"):
                current_pfam_id = line[6:].split(".")[0]
            if line.startswith("//"):
                if current_pfam_id in drift_set:
                    replacement_name = current_hmm_data[1].rstrip()
                    replacement_name += f"|{current_pfam_id}\n"
                    current_hmm_data[1] = replacement_name
                    for output_ln in current_hmm_data:
                        fhOut.write(output_ln)
                current_hmm_data = []
                current_pfam_id = ''
    fhOut.close()


def runhmmemit(num_seqs):
        hmmemit_args = [
            '/home/dbuchan/Applications/hmmer-3.3.2/src/hmmemit',
            '-N',
            num_seqs,
            'hmm_subset.hmm',
        ]
        # print(" ".join(mafft_args))
        try:
            emit_data = subprocess.check_output(hmmemit_args)
            fhOut = open('hmm_generated_seqs.fa', "wb")
            fhOut.write(emit_data)
            fhOut.close()
        except Exception:
            pass


mafft_targets = read_mafft_targets(sys.argv[1])
# print(mafft_targets)
pfam_list = read_summaries(mafft_targets, sys.argv[2])
# print(len(pfam_list))
make_hmms_file(pfam_list, sys.argv[3])
runhmmemit(sys.argv[4])
