import csv
import sys
from collections import defaultdict
import glob
import os
import operator

'''
python calculate_drift_types.py results_data/psiblast_iteration_summaries/
'''

def parse_summary(summary_file):
    family = ''
    families = set()
    drift = False
    erroneous = False
    last_iteration = 0
    data = {}
    summar_reader = None
    with open(summary_file, "r") as fhIn:
        summary_reader = csv.reader(fhIn, delimiter=",")
        try:
            next(summary_reader)
        except:
            erroneous = True
            error_file = summary_file
            return family, drift, erroneous, data
        for row in summary_reader:
            # print(row)
            if int(row[0]) in data:
                data[int(row[0])][row[2]] = int(row[3])
            else:
                data[int(row[0])] = {row[2]: int(row[3])}
            family = row[1]
            families.add(row[2])
            last_iteration = int(row[0])
    if last_iteration != 20:
        erroneous = True
    if len(families) > 1:
        drift = True
    return family, drift, erroneous, data


def calculate_drift_types(main_family, summary):
    print(summary)
    families = set()
    track_data = {}
    for iteration in summary:
        for family in summary[iteration]: 
            families.add(family)
            track_data[family] = {'initial_value': None,
                                  'final_value': None,
                                  'initial_iteration': None,
                                  'final_iteration': None,
                                  'peak_value': None,
                                  }
    for family in families:
        for iteration in summary:
            if family == main_family:
                print(iteration, f"analysing main {family}")
            else:
                print(iteration, f"analysing {family}")

summaries_dir = sys.argv[1]
drift_families = []
non_drift_families = []
erroneous_files = []
for file in glob.glob(f'{summaries_dir}/*.csv'):
    current_family, current_drift, current_erroneous, summary_data = parse_summary(file)
    if current_erroneous:
        erroneous_files.append(file[len(summaries_dir):])
        continue
    if current_drift:
        drift_families.append(current_family)
        calculate_drift_types(current_family, summary_data)
    else:
        non_drift_families.append(current_family)
    break

# print("drifts:", drift_families)
# print("non drifts:", non_drift_families)
# print("errors:", erroneous_files)
# print("drifts:", len(drift_families))
# print("non drifts:", len(non_drift_families))
# print("errors:", len(erroneous_files))
exit()


drift_data = defaultdict(lambda: defaultdict(list))
with open(input_file, "r") as fh:
    driftreader = csv.reader(fh, delimiter=",")
    next(driftreader)
    for row in driftreader:
        # print(row)
        drift_data[row[1]][int(row[2])].append(row[3:])

family_count = 0  # TOTAL NUMBER OF FAMILIES WITH DRIFT
family_hits_sizes = []  # list of the number of drift families in each blast
purified_lost_query = 0
purified_query = 0
spiked_query = 0
grew_query = 0
purified_contaminant = 0
spiked_contaminant = 0
grew_contaminant = 0
purified_total_lost_hit = 0
multiple_families = 0
non_growing_contaminants = 0

for family in drift_data:
    # if family != 'PF05325':
    #     continue
    family_count += 1
    families_set = set()  # The fullset of families seen in this blast run
    final_set = None # set of families in the final iteration
    last_seen_iteration = {} # last iteration we see a family in
    previous_iteration_set = defaultdict(int)  # set of families from the prior
    # iteration
    iteration_set_peaks = defaultdict(int)  # The set of families we've
    # seen in this iteration
    iteration_set_nadirs = defaultdict(int)  # The set of families we've
    # seen in this iteration
    iteration_set_initial = defaultdict(int)  # The set of families we've
    # seen in this iteration
    iteration_set_final = defaultdict(int)  # The set of families we've
    # seen in this iteration
    final_iteration = 0
    for iteration in sorted(drift_data[family]):
        final_iteration = iteration
        iteration_set = defaultdict(int)  # The set of families we've seen in
        # this iteration

        for matches in drift_data[family][iteration]:
            families_set.add(matches[0])
            iteration_set[matches[0]] += int(matches[1])
        # do some min and max tests

        # we now have iteration set which are all the matches seen in a run.
        # loop through the matches again and see when they arrive and go

        contaminant_iteration = 0
        final_set = set()
        for matches in drift_data[family][iteration]:
            hit_family = matches[0]
            num_hits = int(matches[1])
            last_seen_iteration[hit_family] = iteration
            final_set.add(hit_family)
            if iteration == 1:
                iteration_set_peaks[hit_family] = 0
                iteration_set_nadirs[hit_family] = num_hits
            if hit_family not in iteration_set_initial.keys():
                iteration_set_initial[hit_family] = num_hits
                iteration_set_final[hit_family] = num_hits

            iteration_set_final[hit_family] = num_hits

            if num_hits > iteration_set_peaks[hit_family]:
                iteration_set_peaks[hit_family] = num_hits
            if num_hits < iteration_set_nadirs[hit_family]:
                iteration_set_nadirs[hit_family] = num_hits

    
        if family not in iteration_set:
            purified_lost_query += 1
            break
        if (family in iteration_set) and (len(previous_iteration_set.keys())
           < len(iteration_set.keys())):
            purified_total_lost_hit += 0
        previous_iteration_set = iteration_set

    for hit_family in families_set:
        if hit_family not in final_set:
            iteration_set_final[hit_family] = 0

    print("family", family)
    print("Peak", iteration_set_peaks)
    print("Lowest", iteration_set_nadirs)
    print("First value", iteration_set_initial)
    print("Final value", iteration_set_final)
    print("Last seen iter", last_seen_iteration)
    number_families_at_end = 0
    small_non_growing_contaminant = False
    query_grew = False
    query_spiked = False
    query_purified = False
    contaminant_grew = False
    contaminant_spiked = False
    contaminant_purified = False
    max_pairing = max(iteration_set_final.items(), key=operator.itemgetter(1))

    for hit_family in families_set:
        if last_seen_iteration[hit_family] == final_iteration:
            number_families_at_end += 1
        ten_percent_of_peak = iteration_set_peaks[hit_family]*0.2

        if hit_family == family:
            if iteration_set_peaks[hit_family] \
               > iteration_set_initial[hit_family]+ten_percent_of_peak:
                query_grew = True
            if iteration_set_final[hit_family] \
               < iteration_set_peaks[hit_family]*0.8:
                query_spiked = True
            if iteration_set_final[hit_family] \
               < iteration_set_peaks[hit_family]*0.2:
                query_purified = True
            if max_pairing[0] != hit_family and max_pairing[1]*0.9 \
               > iteration_set_final[hit_family]:
                query_purified = True
        else:
            if iteration_set_peaks[hit_family] \
               > iteration_set_initial[hit_family]+ten_percent_of_peak:
                contaminant_grew = True
            if iteration_set_peaks[hit_family] \
               > iteration_set_final[hit_family]*0.8:
                contaminant_spiked = True
            if iteration_set_final[hit_family]*2 \
               < iteration_set_peaks[hit_family]*0.2:
                contaminant_purified = True

    for hit_family in families_set:
        if hit_family != family:
            if iteration_set_initial[hit_family] <= max_pairing[1]*0.1 and iteration_set_final[hit_family] <= max_pairing[1]*0.1:
                small_non_growing_contaminant = True
    
    if small_non_growing_contaminant:
        non_growing_contaminants += 1
    if number_families_at_end > 2:
        multiple_families += 1
    if query_purified:
        purified_query += 1
    if query_spiked:
        spiked_query += 1
    if query_grew:
        grew_query += 1
    if contaminant_purified:
        purified_contaminant += 1
    if contaminant_spiked:
        spiked_contaminant += 1
    if contaminant_grew:
        grew_contaminant += 1
    # exit()
    family_hits_sizes.append(len(families_set))

# print(f"Mean families found: {statistics.mean(family_hits_sizes)}")
# print(f"SD families found: {statistics.stdev(family_hits_sizes)}")
print("---")
print(f"Purifying Selection: LOST QUERY: {purified_lost_query}")
print(f"Purifying Selection: QUERY REDUCED: {purified_query}")
print(f"Purifying Selection: CONTAMINANT REDUCED: {purified_contaminant}")
print(f"Purifying Selection: TOTALLY MISSING HIT: {purified_total_lost_hit}")
print(f"Count with Growing Query: {grew_query}")
print(f"Count with Spiked Query: {spiked_query}")
print(f"Count with Growing Contaminant: {grew_contaminant}")
print(f"Count with Spiked Contaminant: {spiked_contaminant}")
print(f"Count with Non growing contaminants: {non_growing_contaminants}")
print(f"Count with multiple contaminants: {multiple_families}")