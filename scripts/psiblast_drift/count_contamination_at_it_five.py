import csv
import glob

#
# python count_contamination_at_it_five.py /home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries
#
query_purified = "set_where_the_query_was_purified_out.txt"
contam_purified = "set_where_contaminants_are_purified_out.txt"
contam_grew = "set_where_contaminants_grew.txt"
contam_complex = "set_with_complex_contamination_behaviours.txt"
files = [query_purified, contam_purified, contam_grew, contam_complex]

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


drift_summaries = read_summaries("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/")
query_purified = "set_where_the_query_was_purified_out.txt"
contam_purified = "set_where_contaminants_are_purified_out.txt"
contam_grew = "set_where_contaminants_grew.txt"
contam_complex = "set_with_complex_contamination_behaviours.txt"
insig = "insignificant_drifts.txt"
non_drift = "non_drift_list.txt"
summaries_location = "/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_drift_summary"

files = [insig, non_drift, query_purified, contam_purified, contam_grew, contam_complex]

# find out ratio at iteration 5
for file in files:
    print(file)
    pfam_set = read_list(summaries_location, file)
    for pfam in pfam_set:
        print(drift_summaries[pfam][1])