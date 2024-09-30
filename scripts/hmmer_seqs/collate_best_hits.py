import glob
import sys
from collections import defaultdict

#
# python scripts/hmmer_seqs/collate_best_hits.py ./results_data/hmm_closest/ > results_data/hmm_closest/drift_summary.csv
#

collation = {}
for file in glob.glob(f'{sys.argv[1]}/*.best'):
    with open(file, "r", encoding="utf-8") as fhIn:
        target_family = ''
        for line in fhIn:
            entries = line.split(",")
            # print(entries)
            target = entries[0].split("-")
            target_family = target[0]

            if target_family not in collation:
                collation[target_family] = {"nobest_count": 0,
                                            "drift_count": 0,
                                            "correct_count": 0,
                                            "total_count": 0}
            if "None" in entries[1]:
                collation[target_family]["nobest_count"] += 1
            elif target_family in entries[1]:
                collation[target_family]["correct_count"] += 1
            else:
                collation[target_family]["drift_count"] += 1
            collation[target_family]["total_count"] += 1


print("target_family,tot_seqs,correct_count,drift_count,nobest")
for family in collation:
    total = collation[family]['total_count']
    correct = collation[family]['correct_count']
    drift = collation[family]['drift_count']
    nobest = collation[family]['nobest_count']
    print(f'{family},{total},{correct},{drift},{nobest}')

