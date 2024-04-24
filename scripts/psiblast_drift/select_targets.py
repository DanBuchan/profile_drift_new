import random

def select_from_list(file, number):
    target_list = []
    with open(file, "r", encoding="utf-8") as fhIn:
        for line in file:
            pf_family = line.rstrip()
            target_list.append(pf_family)
    return(random.sample(target_list, number))
    

print("family,type")
non_drift_list = "results_data/drift_summary/non_drift_list.txt"
non_drift_set = select_from_list(non_drift_list, 100)
for item in non_drift_set:
    print(f"{item},non_drift")