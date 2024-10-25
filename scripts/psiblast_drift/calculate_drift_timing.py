import csv
import glob

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


drift_summaries = read_summaries("/home/dbuchan/Projects/profile_drift/results_data/drift/pfam_rep_psiblast_iteration_summaries/")
# print(drift_summaries)

# calulate mean contaminants
total_psiblasts = 0
total_contaminants = 0
for file in files:
    pfam_set = read_list(summaries_location, file)
    for pfam in pfam_set:
        if pfam in drift_summaries:
            total_psiblasts += 1
            hit_set = set()
            for iteration in drift_summaries[pfam]:
                for hit_pfam in drift_summaries[pfam][iteration]:
                    # print(hit_pfam)
                    hit_set.add(hit_pfam)
            total_contaminants +=len(hit_set) -1
#     break

print(f"Total runs with contaminants: {total_psiblasts}")
print(f"Total number of contaminants: {total_contaminants}")
print(f"Mean number of contaminants: {total_contaminants/total_psiblasts}")

# Query purified out
pfam_set = read_list(summaries_location, query_purified)
total_families = 0
total_iterations = 0
for pfam in pfam_set:
    if pfam in drift_summaries:
        last_seen_iteration = 1
        total_families += 1
        for iteration in drift_summaries[pfam]:
            hit_set = set()
            for hit_pfam in drift_summaries[pfam][iteration]:
                # print(hit_pfam)
                hit_set.add(hit_pfam)
            if pfam not in hit_set:
                last_seen_iteration = iteration-1
                break
        total_iterations += last_seen_iteration
    
print(f"Total runs for query purified: {total_families}")
print(f"Mean iteration the query last seen: {total_iterations/total_families}")

# contaminant purified out
pfam_set = read_list(summaries_location, contam_purified)
# print(len(pfam_set))
# first find the subset with only 2 for easy calc
sub_set = []
for pfam in pfam_set:
    if pfam in drift_summaries:
        max_contamin = 0
        for iteration in drift_summaries[pfam]:
            hit_set = set()
            for hit_pfam in drift_summaries[pfam][iteration]:
                hit_set.add(hit_pfam)
            # print(hit_set)
            if len(hit_set) > max_contamin:
                max_contamin = len(hit_set)

        if max_contamin == 2:
            sub_set.append(pfam)

total_families = 0
tot_first_seen_iter = 0
tot_last_seen_iter = 0
for pfam in sub_set:
    if pfam in drift_summaries:
        last_seen_contamin = 0
        first_seen_contamin = 0
        total_families += 1
        for iteration in drift_summaries[pfam]:
            hit_set = set()
            for hit_pfam in drift_summaries[pfam][iteration]:
                hit_set.add(hit_pfam)
                if hit_pfam not in pfam and first_seen_contamin == 0:
                    first_seen_contamin = iteration
            if len(hit_set) == 1:
                last_seen_contamin = iteration - 1
        
        tot_first_seen_iter += first_seen_contamin
        tot_last_seen_iter += last_seen_contamin

print(f"Total runs for contaminants purified: {total_families}")
print(f"Mean iteration contaminant first seen: {tot_first_seen_iter/total_families}")
print(f"Mean iteration contaminants removed: {tot_last_seen_iter/total_families}")


total_families = 0
total_sig_iteration = 0
for file in ["set_where_contaminants_grew.txt", "set_with_complex_contamination_behaviours.txt"]:
    pfam_set = read_list(summaries_location, file)
    for pfam in pfam_set:
        if pfam in drift_summaries:
            total_families += 1
            sig_iteration = 0
            for iteration in drift_summaries[pfam]:
                total_query = 0
                total_contam_hits = 0
                
                for hit_pfam in drift_summaries[pfam][iteration]:
                    # print(pfam, hit_pfam)
                    if pfam in hit_pfam:
                        total_query =  drift_summaries[pfam][iteration][hit_pfam]
                    else:
                        total_contam_hits += drift_summaries[pfam][iteration][hit_pfam]
                if total_query == 0:
                    percentage = 100
                else:
                    percentage = (total_contam_hits/total_query) * 100
                if sig_iteration == 0 and percentage >= 5:
                    sig_iteration = iteration
            total_sig_iteration += sig_iteration

print(f"Total runs for complex contaminants: {total_families}")
print(f"Mean iteration of signif contamination: {total_sig_iteration/total_families}")
