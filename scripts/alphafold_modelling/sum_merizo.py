import csv
import sys
import statistics

#
# python sum_merizo.py results_data/alphafold_models/merizo_search_hits.csv
#

def read_merizo(file):
    data = {}
    with open(file, "r", encoding="utf-8") as fh:
        next(fh)
        merizoreader = csv.reader(fh, delimiter=",")
        for row in merizoreader:
            if row[0] not in data:
                data[row[0]] = {}
            if row[1] not in data[row[0]]:
                data[row[0]][row[1]] = {}
            data[row[0]][row[1]][int(row[2])] = {'cath_code': row[5], 'tm_score': row[4]}
            # print(row)
    return(data)

def read_counts(data, iterations):
    count_NAs = 0
    count_matched = 0
    count_matched_unassigned = 0
    
    count_became_unassigned = 0
    count_became_na = 0
    count_missed = 0
    first_iteration_tm_scores = []
    final_iteration_tm_scores = []
    for prot_id in data:
        cath_codes = set()
        if 1 not in data[prot_id] or 20 not in data[prot_id]:
            continue
        for iteration in iterations:
            cath_codes.add(data[prot_id][iteration]['cath_code'])
        if 'NA' not in data[prot_id][1]['tm_score']:
            first_iteration_tm_scores.append(float(data[prot_id][1]['tm_score']))
        if 'NA' not in data[prot_id][20]['tm_score']:
            final_iteration_tm_scores.append(float(data[prot_id][20]['tm_score']))

        if len(cath_codes) == 1:
            if 'NA' in cath_codes:
                count_NAs += 1
            elif 'UNASSIGNED' in cath_codes:
                count_matched_unassigned += 1
            else:
                count_matched += 1
        else:
            if 'NA' in cath_codes:
                count_became_na += 1
            elif 'UNASSIGNED' in cath_codes:
                count_became_unassigned += 1
            else:
                count_missed += 1
    print(f'Total Analysed {len(data)}')
    print(f'No domains recognised: {count_NAs}')
    print(f'Domain always the same: {count_matched}')
    print(f'Domain always unassigned: {count_matched_unassigned}')
    print(f'Domain id changed: {count_missed}')
    print(f'Domain changed to NA: {count_became_na}')
    print(f'Domain changed to unassigned: {count_became_unassigned}')
    print(f'mean tmscore at 1: {statistics.mean(first_iteration_tm_scores)}')
    print(f'mean tmscore at 10: {statistics.mean(final_iteration_tm_scores)}')
    
    


types = ['contaminants_complex', 'contaminants_grew', 'contaminants_purified',
         'insig_drift', 'non_drift', 'query_purified']


merizo_data = read_merizo(sys.argv[1])

c_complex_data = merizo_data['contaminants_complex']
c_grew_data = merizo_data['contaminants_grew']
c_purified_data = merizo_data['contaminants_purified']
insig_data = merizo_data['insig_drift']
non_data = merizo_data['non_drift']
q_purified_data = merizo_data['query_purified']

print("Non Drift")
read_counts(non_data, [1, 20])
print("Insig Drift")
read_counts(insig_data, [1, 20])
print("Query purified")
read_counts(q_purified_data, [1, 20])
print("Contaminants Grew")
read_counts(c_grew_data, [1, 20])
print("Contaminants Purified")
read_counts(c_purified_data, [1, 20])
print("Contaminants Complex")
read_counts(c_complex_data, [1, 20])

#print(merizo_data)