import sys
import csv
from collections import defaultdict
import glob
from subprocess import Popen, PIPE
import os

# This script untars the alignments and preps a fasta file and an a3m file for alphafold to do it's work

# 1. read in target list to get class
# for each target
# 2. unpack the seq tarball
# 3. work out which set of seqs need to be aligned as pre the readme rules
# 4. use clustalOmega to align the seqs outputting in aln or a3m format 

# 1. non_drift - 2 models, 1st and the last : 100 examples == 200 models
# 2. insig_drift - 2 models build the 1st and the last : 100 examples == 200 models
# 3. contaminants_grew - 4 models, 1st, last and 5th and 15th : 100 examples == 400 models
# 4. contaminants_purified - 3 models 1st, last and peak contamination : 50 examples = 150 models
# 5. contaminants_complex -  3 models 1st, last and peak contamination : 100 examples = 300 models *
# 6. query_purified - 4 models, 1st, last and 5th and 15th : 100 examples == 400 models

# might be sensible to arrange as sets given the type and then run multiple array jobs of 100 on the cluster

# usage: python prep_alignments.py /home/dbuchan/Projects/profile_drift/results_data/generation_or_af_targets/alphafold_targets.csv /home/dbuchan/Projects/profile_drift/results_data/alphafold_targets/ 

# /home/dbuchan/Projects/profile_drift/results_data/alphafold_targets/sorted_target_data

def read_target_sample(file):
    with open(file, encoding="utf-8") as fhIn:
        next(fhIn)
        lookup = defaultdict(list)
        targetreader = csv.reader(fhIn, delimiter=',')
        for row in targetreader:
            # print(row)
            if len(row) > 0:
                lookup[row[0]].append(row[1])
    return lookup

def untar_file(file):
    args = ['/usr/bin/tar',
            '-zxvf',
            f'{file}'
    ]
    # print("Untarring", " ".join(args))
    try:
        p = Popen(args, stdout=PIPE, stderr=PIPE)
        result_stdout, err = p.communicate()
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if p.returncode != 0:
        print("Non Zero Exit status: "+str(p.returncode))
        raise OSError("Non Zero Exit status: "+str(p.returncode))

def create_fasta(main_pfam, location):
    for file in glob.glob("*_iteration1_*"):
        # print(file)
        with open(file, "r") as fhIn:
            header = next(fhIn)
            seq = next(fhIn)
        fhOut = open(f'{sys.argv[2]}sorted_target_data/{location}/{main_pfam}.fa', 'w', encoding="utf-8")
        fhOut.write(header)
        fhOut.write(seq)
        fhOut.close()
      
def align_seqs(align_list, location, main_pfam):
    for iteration in align_list:
        # print(iteration)
        for file in glob.glob(f'*_iteration{iteration}_*'):
            args = ['/home/dbuchan/Applications/clustalo_1_2_4/clustalo',
                    '--in',
                    file,
                    '--out',
                    f'{sys.argv[2]}sorted_target_data/{location}/{main_pfam}_{iteration}.a2m',
                    '--outfmt=a2m',
                    '--threads=10',
            ]
            print("Clustering", " ".join(args))
            try:
                p = Popen(args, stdout=PIPE, stderr=PIPE)
                result_stdout, err = p.communicate()
            except Exception as e:
                print(str(e))
                sys.exit(1)
            if p.returncode != 0:
                print("Non Zero Exit status: "+str(p.returncode))
                ## raise OSError("Non Zero Exit status: "+str(p.returncode))

def process_complex(target_data, main_pfam):
    # 3 models 1st, last and peak contamination : 100 examples = 300 models
    peak_iteration = 1
    max_contamination_percentage = 0
    last_iteration = 0
    for iteration in list(target_data.keys()):
        family_count = 0
        contaminant_value = 0
        for family in target_data[iteration]:
            #print(iteration)
            try:
                family_count = int(target_data[iteration][main_pfam])
            except Exception as e:
                family_count = 0
            if family not in main_pfam:
                contaminant_value += int(target_data[iteration][family])
        if family_count == 0:
            contamination_percentage = contaminant_value
        else:
            contamination_percentage = contaminant_value/family_count
        if contamination_percentage > max_contamination_percentage:
            max_contamination_percentage = contamination_percentage
            peak_iteration = iteration
        last_iteration = int(iteration)
        # print(iteration, family_count, contaminant_value, contamination_percentage)
        # print(peak_iteration)
    create_fasta(main_pfam, "contaminants_complex")
    align_seqs([1, peak_iteration, last_iteration], "contaminants_complex", main_pfam)

def process_purified(target_data, main_pfam):
    # 3 models 1st, last and peak contamination : 100 examples = 300 models
    peak_iteration = 1
    max_contamination_percentage = 0
    last_iteration = 0
    for iteration in list(target_data.keys()):
        family_count = 0
        contaminant_value = 0
        for family in target_data[iteration]:
            try:
                family_count = int(target_data[iteration][main_pfam])
            except Exception as e:
                family_count = 0
            if family not in main_pfam:
                contaminant_value += int(target_data[iteration][family])
        if family_count == 0:
            contamination_percentage = contaminant_value
        else:
            contamination_percentage = contaminant_value/family_count
        if contamination_percentage > max_contamination_percentage:
            max_contamination_percentage = contamination_percentage
            peak_iteration = iteration
        last_iteration = int(iteration)
        # print(iteration, family_count, contaminant_value, contamination_percentage)
        # print(peak_iteration)
    create_fasta(main_pfam, "contaminants_purified")
    align_seqs([1, peak_iteration, last_iteration], "contaminants_purified", main_pfam)


def process_non_drift(target_data, main_pfam):
    create_fasta(main_pfam, "non_drift")
    align_seqs([1, 20], "non_drift", main_pfam)

def process_insig(target_data, main_pfam):
    create_fasta(main_pfam, "insig_drift")
    align_seqs([1, 20], "insig_drift", main_pfam)

def process_grew(target_data, main_pfam):
    create_fasta(main_pfam, "contaminants_grew")
    align_seqs([1, 5, 15, 20], "contaminants_grew", main_pfam)

def process_q_purified(target_data, main_pfam):
    create_fasta(main_pfam, "query_purified")
    align_seqs([1, 5, 15, 20], "query_purified", main_pfam)

def tidy_up():
    for file in glob.glob("*.fa"):
        os.remove(file)

# generate the fasta sequence and alignment
def process_targets(af_dir, target_class):
    for file in glob.glob(f"{af_dir}/*.csv"):
        print(file)
        name_stub = file[:-17]
        tarball = f"{name_stub}seqs.tar.gz"
        untar_file(tarball)
        target_family = ''
        target_data = defaultdict(dict)
        with open(file, encoding="utf-8") as fhIn:
            next(fhIn)
            summaryreader = csv.reader(fhIn, delimiter=',')
            for row in summaryreader:
                # print(row)
                target_family = row[1]
                target_data[row[0]][row[2]] = row[3]
        #print(target_family)
        # print(target_data)
        
        for target_type in target_class[target_family]:
            print(target_type)
            # if "contaminants_complex" in target_type:
            #     continue
            if "contaminants_complex" in target_type:
                process_complex(target_data, target_family)
            if "non_drift" in target_type:
                process_non_drift(target_data, target_family)
            if "insig_drift" in target_type:
                process_insig(target_data, target_family)
            if "contaminants_grew" in target_type:
                process_grew(target_data, target_family)
            if "contaminants_purified" in target_type:
                process_purified(target_data, target_family)
            if "query_purified" in target_type:
                process_q_purified(target_data, target_family)
                
        tidy_up()
        # break


target_class = read_target_sample(sys.argv[1])
print(target_class)
process_targets(sys.argv[2], target_class)



