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
    
    # print(track_data)
    for family in track_data.keys():
        if family in main_family:
            continue # don't compare the main family
        
        try: # if a search never finds it's own family set this threshold to 0
            five_percent_peak = track_data[main_family]['peak_value'] * 0.05
        except:
            five_percent_peak = 0

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
    # if 'results_data/psiblast_iteration_summaries/10047_GH123_N-consensus_blast_summary.csv' not in file:
    #     continue
    # print(file)
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
    # if count == 2:
    #     break

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
contaminants_purified  = set()
contaminants_grew  = set()
contaminants_spiked  = set()
families_lost  = set()
contaminants_stable  = set()
contaminants_complex  = set()
query_purified = []
query_grew = []
query_spiked = []
query_stable = []


for family in significant_drifts:
    count+=1
    data = full_results[family]
    # print(family)
    # pprint.pp(data)
    if data['families_at_final_iteration'] < data['number_of_drift_familes'] + 1:
        families_lost.add(family)
    growth_data = data['growth_types']
    grow_set = set()
    spiked_set = set()
    purified_set = set()
    flat_set = set()
    # pprint.pp(growth_data)
    # First total up what the query family is doing
    if family in growth_data:
        if growth_data[family]['grew']:
            query_grew.append(family)
        if growth_data[family]['spiked']:
            query_spiked.append(family)
        if growth_data[family]['flat']:
            query_stable.append(family)
        if growth_data[family]['purified']:
            query_purified.append(family)
        del growth_data[family]
    
    # loop over the contaminants and check what the non-neligible ones are doing.
    # Add the behaviours to theses sets
    for contaminant_family in growth_data:
        if growth_data[contaminant_family]['negligible_contaminant'] is False:
            grow_set.add(growth_data[contaminant_family]['grew'])
            spiked_set.add(growth_data[contaminant_family]['spiked'])
            purified_set.add(growth_data[contaminant_family]['purified'])
            flat_set.add(growth_data[contaminant_family]['flat'])

    if len(grow_set) == 1 and len(spiked_set) == 1 and len(purified_set) == 1 and len(flat_set) == 1: 
        if list(spiked_set)[0] or list(grow_set)[0]:
            contaminants_grew.add(family)
        if list(spiked_set)[0]:
            contaminants_spiked.add(family)
        if list(flat_set)[0]:
            contaminants_stable.add(family)
        if list(purified_set)[0]:
            contaminants_purified.add(family)
        if (list(spiked_set)[0] or list(grow_set)[0]) and list(flat_set)[0]:
            contaminants_complex.add(family)
        if (list(spiked_set)[0] or list(grow_set)[0]) and list(purified_set)[0]:
            contaminants_complex.add(family)
            
    else:
        contaminants_complex.add(family)

print(f"Total Analysed drift familes: {count}")
print(f"Count where there were fewer families at iteration 20 than total seen: {len(list(families_lost))}")
print(f"Count where query grew: {len(query_grew)}")
print(f"Count where query peaked: {len(query_spiked)}")
print(f"Count where query had flat growth: {len(query_stable)}")
print(f"Count where query was purified: {len(query_purified)}")
print(f"Count where contaminants with different behaviours: {len(list(contaminants_complex))}")
print(f"Count where contaminants only grew: {len(list(contaminants_grew))}")
print(f"Count where contaminants peaked: {len(list(contaminants_spiked))}")
print(f"Count where contaminants only had flat growth: {len(list(contaminants_stable))}")
print(f"Count where contaminants only purified: {len(list(contaminants_purified))}")

with open("set_where_there_are_fewer_families_at_the_end.txt", "w", encoding="utf-8") as fhOut:
    for item in families_lost:
        fhOut.write(f'{item}\n')

with open("set_where_the_query_was_purified_out.txt", "w", encoding="utf-8") as fhOut:
    for item in query_purified:
        fhOut.write(f'{item}\n')

with open("set_with_complex_contamination_behaviours.txt", "w", encoding="utf-8") as fhOut:
    for item in contaminants_complex:
        fhOut.write(f'{item}\n')

with open("set_where_contaminants_are_purified_out.txt", "w", encoding="utf-8") as fhOut:
    for item in contaminants_purified:
        fhOut.write(f'{item}\n')

with open("set_where_contaminants_are_purified_out.txt", "w", encoding="utf-8") as fhOut:
    for item in contaminants_purified:
        fhOut.write(f'{item}\n')

with open("set_where_contaminants_grew.txt", "w", encoding="utf-8") as fhOut:
    for item in contaminants_grew:
        fhOut.write(f'{item}\n')