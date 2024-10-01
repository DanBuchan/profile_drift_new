import glob
import sys

#
# python scripts/alphafold_modelling/sum_plddts.py ./results_data/alphafold_models/ > ./results_data/alphafold_models/plddt_summary.csv
#

results = {}
for file in glob.glob(f'{sys.argv[1]}/*.pdb'):
    with open(file, "r", encoding="utf-8") as fh:
        if 'contaminants_complex_PF00106_1_428ec_unrelaxed_rank_1_model_3.pdb' not in file:
            continue
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
                print(f'{res_id}')
                plDDT = line[61:66]
                plDDT = float(plDDT.replace(" ", ""))
                print(f'{plDDT}')
                plDDT_scores[res_id] = plDDT
        print(plDDT_scores)
        plDDT_tot = 0
        plDDT_count = 0
        for res_id in plDDT_scores:
            plDDT_count += 1
            plDDT_tot += plDDT_scores[res_id]
        ave_plDDT = plDDT_tot/res_id
        if drift_class not in results:
            results[drift_class] = {}
        if family not in results[drift_class]:
            results[drift_class][family] = {}
        
        results[drift_class][family][iteration] = ave_plDDT
        print(results)
    # break

print("drift_class,family,iteraton,mean_plDDT")
for drift_class in results:
    for family in results[drift_class]: 
        for iteration in results[drift_class][family]:
            plDDT = results[drift_class][family][iteration]
            print(f'{drift_class},{family},{iteration},{plDDT}')