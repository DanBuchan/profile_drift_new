import glob
import sys

#
# python scripts/alphafold_modelling/sum_plddts.py ./results_data/alphafold_models/*.pdb
#

results = {}
for file in glob.glob(sys.argv[1]):
    with open(file, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startwith("ATOM"):
                entries = line.split()
                print(entries)
    break