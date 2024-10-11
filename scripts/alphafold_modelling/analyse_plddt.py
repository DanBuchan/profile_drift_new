import csv
import sys
#
# python scripts/alphafold_modelling/analyse_plddt.py ./results_data/alphafold_models/plddt_summary.csv  > ./results_data/alphafold_models/plddt_means.csv
#

results = {'contaminants_grew': {},
           'contaminants_complex': {},
           'contaminants_purified': {},
           'insig_drift': {},
           'non_drift': {},
           'query_purified': {},
}
classes = set()
with open(sys.argv[1], "r", encoding="utf-8") as fhIn:
    plddtreader = csv.reader(fhIn, delimiter=',')
    next(plddtreader)
    for row in plddtreader:
        classes.add(row[0])
        if int(row[2]) in results[row[0]]:
            results[row[0]][int(row[2])]["total_plddt"]+=float(row[3])
            results[row[0]][int(row[2])]["class_count"]+=1
        else:
            results[row[0]][int(row[2])] = {"total_plddt": float(row[3]),
                                            "class_count": 1}

for drift_class in results:
    # print(drift_class)
    # print(results[drift_class])
    if "contaminants_grew" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_5 = results[drift_class][5]['total_plddt']/results[drift_class][5]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        print(f'{drift_class},{mean_1},{mean_5},{mean_20}')
    if "query_purified" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_5 = results[drift_class][5]['total_plddt']/results[drift_class][5]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        print(f'{drift_class},{mean_1},{mean_5},{mean_20}')
    if "insig_drift" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        print(f'{drift_class},{mean_1},na,{mean_20}')
    if "non_drift" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        print(f'{drift_class},{mean_1},na,{mean_20}')
    if "contaminants_complex" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        tot_plddt = 0
        tot_count = 0
        for iteration in results[drift_class]:
            if iteration == 1 or iteration == 20:
                continue
            tot_plddt += results[drift_class][iteration]['total_plddt']
            tot_count +=  results[drift_class][iteration]['class_count']
        mean_peak = tot_plddt/tot_count
        print(f'{drift_class},{mean_1},{mean_peak},{mean_20}')
    if "contaminants_purified" in drift_class:
        mean_1 = results[drift_class][1]['total_plddt']/results[drift_class][1]['class_count']
        mean_20 = results[drift_class][20]['total_plddt']/results[drift_class][20]['class_count']
        tot_plddt = 0
        tot_count = 0
        for iteration in results[drift_class]:
            if iteration == 1 or iteration == 20:
                continue
            tot_plddt += results[drift_class][iteration]['total_plddt']
            tot_count +=  results[drift_class][iteration]['class_count']
        mean_peak = tot_plddt/tot_count
        print(f'{drift_class},{mean_1},{mean_peak},{mean_20}')
    