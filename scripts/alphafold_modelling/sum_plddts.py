import glob
import sys

#
# python scripts/alphafold_modelling/sum_plddts.py ./results_data/alphafold_models/
#

results = {}
for file in glob.glob(f'{sys.argv[1]}/*.pdb'):
    with open(file, "r", encoding="utf-8") as fh:
        file_parts = file.split("_")
        print(file_parts)
        drift_class = f'{file_parts[0]}_{file_parts[1]}'
        family = file_parts[2]
        iteration = file_parts[3]

        for line in fh:
            if line.startswith("ATOM"):
                entries = line.split()
                # print(entries)
    break