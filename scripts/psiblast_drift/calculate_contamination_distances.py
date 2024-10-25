import csv
import glob
import numpy as np

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

counts, bins = np.histogram(dist_matrix[dist_matrix>0], bins=np.arange(426463801))
mode_value = np.argmax(counts)

print(f'mode: {mode_value}, mean: {dist_matrix.mean()}, std: {dist_matrix.std()}')