import sys
import numpy as np
import glob
import csv
import pickle

### Usage
# python ./scripts/rep_distance_matrix/build_rep_distance_matrix.py /home/dbuchan/Data/pfam/ ./results_data/distance_matrix/pfam_rep_all_against_all/
###

# open the first file in the dir to get the ordered list of IDs
current_id = ''
start = True
dom_list = []
for file in sorted(glob.glob(f'{sys.argv[1]}*_consensus')):
    # print(file[len(sys.argv[1]):])
    with open(file, "r", encoding="utf-8") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                line = line[1:]
                line = line.rstrip()
                dom_list.append(line)

print(len(dom_list))

similarity_matrix = np.empty((len(dom_list), len(dom_list)), dtype=float)

# populate with values from file
for file in sorted(glob.glob(f'{sys.argv[2]}*.csv')):
    print(file[len(sys.argv[2]):])
    with open(file, 'r') as pffile:
        reader = csv.reader(pffile, delimiter=',')
        next(reader)
        for row in reader:
            try:
                x = dom_list.index(row[0])-1
                y = dom_list.index(row[1])-1
            except Exception:
                continue
            # print(x, y)
            similarity_matrix[x][y] = row[2] 

# now scale/normalise to 0 and 1, flip and fill diagonal
np.fill_diagonal(similarity_matrix, (np.max(similarity_matrix)*1.2))
normalized_similarity_matrix = (similarity_matrix-np.min(similarity_matrix))/(np.max(similarity_matrix)-np.min(similarity_matrix))
flipped_sim_matrix = 1.0 - normalized_similarity_matrix
# now save out distance matrix and dom_list
with open("hmm_rep_distance_matrix.npy", "wb") as f:
    np.save(f, flipped_sim_matrix)
    np.save(f, dom_list)