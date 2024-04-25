import sys

def read_target_id_list(file):
    targets = []
    with open(file, "r", encoding="utf-8") as fhIn:
        for line in fhIn:
            targets.append(line.rstrip())
    return(targets)

target_id = int(sys.argv[1]) # id in the target list to process
targets_file = sys.argv[2] # target list

targets = read_target_id_list(targets_file)

print(targets[target_id])