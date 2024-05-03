import sys
import numpy as np
import pandas as pd
import glob
import csv

### Usage
# python ./scripts/rep_distance_matrix/build_rep_distance_matrix.py ./results_data/distance_matrix/pfam_rep_all_against_all/
###

# open the first file in the dir to get the ordered list of IDs
current_id = ''
start = True
id_list = []
with open(f'{sys.argv[1]}1_hits.csv', "r", encoding="utf-8") as fhIn:
    reader = csv.reader(fhIn, delimiter=',')
    for row in reader:
        if start:
            current_id = row[0]
            id_list.append(row[0])
            start = False
        if row[0] not in current_id:
            print(len(id_list))
            # print(id_list)
            exit()
        id_list.append(row[1])

        
dom_list = set()
# prep the unique list of domains
for file in sorted(glob.glob(f'{sys.argv[1]}*.csv')):
    print(file[len(sys.argv[1]):])
    with open(file, 'r') as pffile:
        reader = csv.reader(pffile, delimiter=',')
        
        for row in reader:
            dom_list.add(row[1])
    break


# df = pd.DataFrame(columns=list(dom_list), index=list(dom_list))
# if 'Imm21-consensus|PF1558' in df.indexes:
#     print("found one")
# if 'ADAMTS_CR_2-consensus|PF17771' in df.columns:
#     print("found two")
# exit()

# for file in sorted(glob.glob(f'{sys.argv[1]}*.csv')):
#     # print(file[len(sys.argv[1]):])
#     with open(file, 'r') as pffile:
#         reader = csv.reader(pffile, delimiter=',')
#         for row in reader:
#             try:
#                 df.loc[row[0], row[1]] = row[2]
#             except Exception as e:
#                 df.loc[row[0], row[1]] = None
#     break

# print(df.max())
# do we need to populate the diagonal?
# for dom in dom_list:
#     df.loc[dom, dom] = 1
 
# print(df)