import glob
import sys

#
# python scripts/alphafold_modelling/sum_plddts.py ./results_data/alphafold_models/
#

results = {}
for file in glob.glob(f'{sys.argv[1]}/*.pdb'):
    with open(file, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                entries = line.split()
                print(entries)
    break