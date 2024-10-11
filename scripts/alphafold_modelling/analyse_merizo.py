import csv
import sys

results = {'contaminants_grew': {},
           'contaminants_complex': {},
           'contaminants_purified': {},
           'insig_drift': {},
           'non_drift': {},
           'query_purified': {},
}
with open(sys.argv[1], "r", encoding="utf-8") as fhIn:
    plddtreader = csv.reader(fhIn, delimiter=',')
    next(plddtreader)
    for row in plddtreader:
        # print(row)
        if row[1] in results[row[0]]:
             results[row[0]][row[1]][int(row[2])]= {"domain_id": row[3],
                                                     "tm_score": row[4],
                                                     "cath_code": row[5]}
        else:
            results[row[0]][row[1]] = {int(row[2]): {"domain_id": row[3],
                                                     "tm_score": row[4],
                                                     "cath_code": row[5]}}
   
print("drift_class,assigned_domains,changed_domain,changed_h_family")
for drift_class in results:
    #print(drift_class)
    changed_domain_count = 0
    changed_cath_count = 0
    domain_count = 0
    for domain in results[drift_class]:
        max_tm = 0
        max_tm_iteration = 0
        min_tm = 100000
        min_tm_iteration = 0
        domain_set = set()
        cath_code_set = set()
        set_count = 0
        analysed = False
        for iteration in results[drift_class][domain]:
            if 'NA' in results[drift_class][domain][iteration]["domain_id"]:
                continue
            analysed = True
            set_count += 1
            domain_set.add(results[drift_class][domain][iteration]["domain_id"])
            cath_code_set.add(results[drift_class][domain][iteration]["cath_code"])
            tm_score = float(results[drift_class][domain][iteration]["tm_score"])
            if tm_score > max_tm:
                max_tm = tm_score
                max_tm_iteration = iteration
            if tm_score < min_tm:
                min_tm = tm_score
                min_tm_iteration = iteration
        if analysed:
            domain_count += 1
        if len(domain_set) > 1:
            changed_domain_count += 1
        if len(cath_code_set) > 1:
            changed_cath_count += 1
    
    print(f'{drift_class},{domain_count},{changed_domain_count},{changed_cath_count}')