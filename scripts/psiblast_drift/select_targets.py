import random

def select_from_list(file, number):
    target_list = []
    with open(file, "r", encoding="utf-8") as fhIn:
        for line in fhIn:
            pf_family = line.rstrip()
            # print(pf_family)
            target_list.append(pf_family)
    if len(target_list) < number:
        # print("hi")
        return target_list
    return random.sample(target_list, number)
    

print("family,type")
drift_list = "results_data/drift_summary/non_drift_list.txt"
drift_set = select_from_list(drift_list, 100)
# print(len(drift_set))
for item in drift_set:
    print(f"{item},non_drift")

drift_list = "results_data/drift_summary/non_drift_list.txt"
drift_set = select_from_list(drift_list, 100)
for item in drift_set:
    print(f"{item},insig_drift")

drift_list = "results_data/drift_summary/set_where_contaminants_grew.txt"
drift_set = select_from_list(drift_list, 100)
for item in drift_set:
    print(f"{item},contaminants_grew")

drift_list = "results_data/drift_summary/set_where_contaminants_are_purified_out.txt"
drift_set = select_from_list(drift_list, 100)
for item in drift_set:
    print(f"{item},contaminants_purified")

# drift_list = "results_data/drift_summary/set_with_complex_contamination_behaviours.txt"
# drift_set = select_from_list(drift_list, 100)
# for item in drift_set:
#     print(f"{item},contaminants_complex")

# drift_list = "results_data/drift_summary/set_where_the_query_was_purified_out.txt"
# drift_set = select_from_list(drift_list, 100)
# for item in drift_set:
#     print(f"{item},query_purified")