import sys
import os
import glob

def read_target_id_list(file):
    targets = []
    with open(file, "r", encoding="utf-8") as fhIn:
        for line in fhIn:
            targets.append(line.rstrip())
    return(targets)

def untar_seqs():
    for tar_file in glob.glob('*.tag.gz')
        tar_args = [
            '/usr/bin/tar',
            'zxf',
            tar_file,
        ]
        tar_output = subprocess.check_output(tar_args)

def align_seqs(seq_file, iterations):
    clean_up = False
    for i, fa_file in enumerate(glob.glob('*.fa')):
        msa = f"{seq_file}_iteration{i+1}_seqs.msa"
        mafft_args = [
        #    '/home/dbuchan/Applications/mafft-7.490-with-extensions/core/mafft',
            '/home/ucbcdwb/Applications/mafft-linux64/mafft.bat',
            fa_file,
        ]
        # print(" ".join(mafft_args))
        try:
            msa_data = subprocess.check_output(mafft_args)
            fhOut = open(msa, "wb")
            fhOut.write(msa_data)
            fhOut.close()
            clean_up=True
        except Exception as e:
            break


target_id = int(sys.argv[1]) # id in the target list to process
targets_file = sys.argv[2] # target list

targets = read_target_id_list(targets_file)
print(f'ALIGNING: {targets[target_id]}')
os.chdir(f'/home/ucbcdwb/Scratch/output/profile/{targets[target_id])}')
untar_seqs()
align_seqs()    
