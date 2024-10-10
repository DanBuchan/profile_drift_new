import glob
import sys
from collections import defaultdict

#
# python scripts/prottrans_seqs/collate_best_hits.py ./results_data/drift/best_hits/prottrans_closest/ > results_data/drift/best_hits/prottrans_closest/drift_summary.csv 
#

collation = {}
for file in glob.glob(f'{sys.argv[1]}/*.best'):
    with open(file, "r", encoding="utf-8") as fhIn:
        target_family = ''
        for line in fhIn:
            entries = line.split(",")
            # print(entries)
            target = entries[0].split("_")
            target_family = target[0]
            mask_level = target[1]

            if target_family not in collation:
                collation[target_family] = {}
            
            if mask_level not in collation[target_family]:
                collation[target_family][mask_level] = {"nobest_count": 0,
                                                        "drift_count": 0,
                                                        "correct_count": 0,
                                                         "total_count": 0}
            if "None" in entries[1]:
                collation[target_family][mask_level]["nobest_count"] += 1
            elif target_family in entries[1]:
                collation[target_family][mask_level]["correct_count"] += 1
            else:
                collation[target_family][mask_level]["drift_count"] += 1
            collation[target_family][mask_level]["total_count"] += 1


print("target_family,mask,tot_seqs,correct_count,drift_count,nobest")
for family in collation:
    for mask in collation[family]:
        total = collation[family][mask]['total_count']
        correct = collation[family][mask]['correct_count']
        drift = collation[family][mask]['drift_count']
        nobest = collation[family][mask]['nobest_count']
        print(f'{family},{mask},{total},{correct},{drift},{nobest}')

