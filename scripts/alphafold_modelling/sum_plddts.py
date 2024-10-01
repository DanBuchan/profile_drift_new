import glob
import sys

#
# python scripts/alphafold_modelling/sum_plddts.py ./results_data/alphafold_models/
#

results = {}
for file in glob.glob(f'{sys.argv[1]}/*.pdb'):
    with open(file, "r", encoding="utf-8") as fh:
        print(file)
        file_parts = file[len(sys.argv[1]):].split("_")
        drift_class = f'{file_parts[0]}_{file_parts[1]}'
        family = file_parts[2]
        iteration = file_parts[3]
        print(drift_class, family, iteration)
        plDDT_scores = {}
        for line in fh:
            if line.startswith("ATOM"):
                res_id = line[24:28]
                res_id = int(res_id.replace(" ", ""))
                # print(f'{res_id}_')
                plDDT = line[61:66]
                plDDT = float(plDDT.replace(" ", ""))
                # print(f'{plDDT}')
                plDDT_scores[res_id] = plDDT
        plDDT_tot = 0
        plDDT_count = 0
        for res_id in plDDT_scores:
            plDDT_count += 1
            plDDT_tot += plDDT_scores[res_id]
        ave_plDDT = plDDT_tot/res_id
        print(ave_plDDT)
    break