import csv
import sys
from collections import defaultdict
import glob
import pprint
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

def return_growth_types(initial, peak, final):
    ten_percent_initial = initial * 0.1
    twenty_percent_initial = initial * 0.2
    eighty_percent_peak = peak * 0.8

    grew = False
    spiked = False
    purified = False
    flat = False

    if peak >= initial + twenty_percent_initial:
        grew = True

    if final <= eighty_percent_peak:
        spiked = True
        grew = True
    
    if final <= initial + ten_percent_initial and grew:
        purified = True

    if grew == False and spiked == False and purified == False:
        flat = True
    
    # and just return flat if abs values are small
    if initial < 30 and peak < 30 and final < 30:
        grew = False
        spiked = False
        purified = False
        flat = True

    return [grew, spiked, purified, flat]


def calculate_drift_types(main_family, summary):
    # summary[iteration][family][value]
    # pprint.pp(summary)
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
                                  'peak_iteration': None
                                  }
    number_families_at_final_iteration = 0
    for family in families:
        for iteration in summary:
            if family in summary[iteration]:
                number_families_at_final_iteration = len(summary[iteration])
                if track_data[family]['initial_iteration'] is None:
                    track_data[family]['initial_iteration'] = iteration
                    track_data[family]['initial_value'] = summary[iteration][family]
                    track_data[family]['peak_value'] = summary[iteration][family]
                    track_data[family]['peak_iteration'] = iteration
                    track_data[family]['final_iteration'] = iteration
                    track_data[family]['final_value'] = summary[iteration][family]
                else:
                    track_data[family]['final_value'] = summary[iteration][family]
                    track_data[family]['final_iteration'] = iteration
                    if summary[iteration][family] > track_data[family]['peak_value']:
                        track_data[family]['peak_value'] = summary[iteration][family]
                        track_data[family]['peak_iteration'] = iteration
                                  
    # pprint.pp(track_data)
    # Now work out drift types
    number_of_contaminents = len(track_data)-1
    results = {'number_of_drift_familes': number_of_contaminents,
               'families_at_final_iteration': number_families_at_final_iteration,
               'growth_types': {},
               }
    for family in track_data.keys():
        drift_types = return_growth_types(track_data[family]['initial_value'],
                                          track_data[family]['peak_value'],
                                          track_data[family]['final_value'])
        results['growth_types'][family] = {'grew': drift_types[0],
                                           'spiked': drift_types[1],
                                           'purified': drift_types[2],
                                           'flat': drift_types[3]}
        
    for family in track_data.keys():
        if family in main_family:
            continue # don't compare the main family
        five_percent_peak = track_data[main_family]['peak_value'] * 0.05
        if track_data[family]['peak_value'] <= five_percent_peak and \
               track_data[family]['final_value'] <= five_percent_peak:
           results['growth_types'][family]['negligible_contaminant'] = True
        else:
            results['growth_types'][family]['negligible_contaminant'] = False
    return(results)
   

summaries_dir = sys.argv[1]
drift_families = []
non_drift_families = []
erroneous_files = []
full_results = {}
count = 0
for file in glob.glob(f'{summaries_dir}/*.csv'):
    count += 1
    current_family, current_drift, current_erroneous, summary_data = parse_summary(file)
    if current_erroneous:
        erroneous_files.append(file[len(summaries_dir):])
        continue
    if current_drift:
        drift_families.append(current_family)
        full_results[current_family] = calculate_drift_types(current_family, summary_data)
    else:
        non_drift_families.append(current_family)
    if count == 40:
        break

###
### Output lists of results
###
        
with open("drift_list.txt", "w") as fhDrifts:
    for family in drift_families:
        fhDrifts.write(f'{family}\n')

with open("non_drift_list.txt", "w") as fhNDrifts:
    for family in non_drift_families:
        fhNDrifts.write(f'{family}\n')

with open("drift_error_list.txt", "w") as fhErrors:
    for family in erroneous_files:
        fhErrors.write(f'{family}\n')

negligible_drifts = []
significant_drifts = []
significant_drift = None
for main_family in full_results.keys():
    significant_drift = False
    drift_info = full_results[main_family]['growth_types']
    for family in drift_info:
        if family in main_family:
            continue
        if full_results[main_family]['growth_types'][family]['negligible_contaminant'] is False:
            significant_drift = True
    if significant_drift:
        significant_drifts.append(main_family)
    else:
        negligible_drifts.append(main_family)

with open("significant_drifts.txt", "w") as fhSigDrifts:
    for family in significant_drifts:
        fhSigDrifts.write(f'{family}\n')

with open("insignificant_drifts.txt", "w") as fhInSigDrifts:
    for family in negligible_drifts:
        fhInSigDrifts.write(f'{family}\n')

###
### Output summary
###
count = 0
contaminants_purified = 0
contaminants_grew = 0
contaminants_spiked = 0
families_lost = 0
contaminants_stable = 0
contaminants_complex = 0
query_purified = 0
query_grew = 0
query_spiked = 0
query_lost = 0
query_stable = 0

for family in significant_drifts:
    count+=1
    data = full_results[family]
    print(family)
    pprint.pp(data)
    if data['families_at_final_iteration'] < data['number_of_drift_familes'] + 1:
        families_lost += 1
    growth_data = data['growth_types']
    grow_set = set()
    spiked_set = set()
    purified_set = set()
    flat_set = set()
    for growth_family in growth_data:
        if family in growth_family:
            if growth_data[family]['grew']:
                query_grew += 1
            if growth_data[family]['spiked']:
                query_spiked += 1
            if growth_data[family]['flat']:
                query_stable+= 1
            if growth_data[family]['purified']:
                query_purified+= 1
        elif growth_data[growth_family]['negligible_contaminant'] == False:
            grow_set.add(growth_data[family]['grew'])
            spiked_set.add(growth_data[family]['grew'])
            purified_set.add(growth_data[family]['grew'])
            flat_set.add(growth_data[family]['grew'])
    if len(grow_set) == 1 and len(spiked_set) == 1 and len(purified_set) == 1 and len(flat_set) == 1: 
        if list(spiked_set)[0] or list(grow_set)[0]:
            contaminants_grew +=1
        if list(spiked_set)[0]:
            contaminants_spiked+=1
        if list(flat_set)[0]:
            contaminants_stable+=1
        if list(purified_set)[0]:
            contaminants_purified+=1
    else:
        contaminants_complex+=1

print(f"Total Analysed familes: {count}")
print(f"Count where there were fewer families at iteration 20 than total seen: {families_lost}")
print(f"Count where query grew: {query_grew}")
print(f"Count where query peaked: {query_spiked}")
print(f"Count where query had flat growth: {query_stable}")
print(f"Count where query was purified: {query_purified}")
print(f"Count where contaminants only grew: {contaminants_complex}")
print(f"Count where contaminants peaked: {contaminants_complex}")
print(f"Count where contaminants only stable: {contaminants_complex}")
print(f"Count where contaminants only purified: {contaminants_complex}")


# print(f"Purifying Selection: LOST QUERY: {purified_lost_query}")
# print(f"Purifying Selection: QUERY REDUCED: {purified_query}")
# print(f"Purifying Selection: CONTAMINANT REDUCED: {purified_contaminant}")
# print(f"Purifying Selection: TOTALLY MISSING HIT: {purified_total_lost_hit}")
# print(f"Count with Growing Query: {grew_query}")
# print(f"Count with Spiked Query: {spiked_query}")
# print(f"Count with Growing Contaminant: {grew_contaminant}")
# print(f"Count with Spiked Contaminant: {spiked_contaminant}")
# print(f"Count with Non growing contaminants: {non_growing_contaminants}")
# print(f"Count with multiple contaminants: {multiple_families}")