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
        for line in fh:
            if line.startswith("ATOM"):
                plDDT = line[61:66]
                print(f'{plDDT}_')
    break