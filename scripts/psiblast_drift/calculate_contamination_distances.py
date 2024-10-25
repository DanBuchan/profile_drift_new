import csv
import glob
import numpy as np
# from matplotlib import pyplot as plt

def read_summaries(path):
    summaries = {}
    for file in glob.glob(f'{path}*'):
        # if "19377_A0A8C4KJU7.1_9-26_blast_summary.csv" not in file:
        #     continue
        # print(file)
        with open(file, "r", encoding="utf-8") as fh:
            next(fh)
            summaryreader = csv.reader(fh, delimiter=',')
            for i, row in enumerate(summaryreader):
                if i == 0:
                    summaries[row[1]] = {}
                if int(row[0]) not in summaries[row[1]]:
                    # print(row)
                    summaries[row[1]][int(row[0])] = {}
                # print(summaries)
                summaries[row[1]][int(row[0])][row[2]] = int(row[3])
    return(summaries)

def read_list(path, file):
    pfam_list = []
    with open(f'{path}/{file}', "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            pfam_list.append(line)
    return(pfam_list)

summaries_location = "/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_drift_summary"
# non_drift = "non_drift_list.txt"
# insig_drift = "insignificant_drifts.txt"
query_purified = "set_where_the_query_was_purified_out.txt"
contam_purified = "set_where_contaminants_are_purified_out.txt"
contam_grew = "set_where_contaminants_grew.txt"
contam_complex = "set_with_complex_contamination_behaviours.txt"
files = [query_purified, contam_purified, contam_grew, contam_complex]

with open('/home/dbuchan/Projects/profile_drift/results_data/distance_matrix/pfam_rep_all_against_all/rand_rep_distance_matrix.npy', 'rb') as f:
    dist_matrix = np.load(f)
    dom_list = np.load(f)
drift_summaries = read_summaries("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/")

# bins = list(np.arange(0.0,1.0,0.001))
# flat_array = dist_matrix.flatten()
# for value in flat_array:
#     print(value)
# frq, edges = np.histogram(flat_array, bins)
# fig, ax = plt.subplots()
# ax.bar(edges[:-1], frq, width=np.diff(edges), edgecolor="black", align="edge")

# plt.show()

# print(f'mean: {dist_matrix.mean()}, std: {dist_matrix.std()}')

total_dist = 0
total_comparisions =0
for file in files:
    pfam_set = read_list(summaries_location, file)
    pfam_id_set = set()
    for pfam in pfam_set:
        pfam_id_set = set()
        if pfam in drift_summaries:
            for iteration in drift_summaries[pfam]:
                for hit_pfam in drift_summaries[pfam][iteration]:
                    if pfam not in hit_pfam:
                        pfam_id_set.add(hit_pfam)
        query_idx = dom_list.index(pfam)
        for hit_pfam in pfam_id_set:
            hit_idx = dom_list.index(hit_pfam)
            print(dist_matrix[query_idx,hit_idx])
            exit