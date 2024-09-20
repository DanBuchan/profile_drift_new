#!/bin/bash -l
#$ -o /home/ucbcdwb/Scratch/output/alpha/std.out
#$ -e /home/ucbcdwb/Scratch/output/alpha/std.err

# Request an hour of run time
#$ -l h_rt=200:00:0

# Request 10 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=10G

#$ -N runalpha
#$ -l tmem=2G
#$ -l gpu=true
#$ -pe gpu 1
#$ -R y

# Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
#$ -t 1-2
# Set the name of the job.

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/cuda-11.8/lib64/
conda activate /cluster/project1/ProCovar/lmoffat/miniconda/envs/alphafold 
cd /cluster/project1/ProCovar/colabfold-customtemplates

CLASS='contaminants_grew' # remove 57 leading chars
FASTA_ID="$(tail -n +${SGE_TASK_ID} /home/dbuchan/inputs/alphafold/${CLASS}/targets | head -n 1)"
FAMILY_ID=${FASTA_ID::-3}

for i in /home/dbuchan/inputs/alphafold/${CLASS}/${FAMILY_ID}_*.a2m; do
    # bash /cluster/project1/ProCovar/colabfold-customtemplates/main_db.sh ${FASTA_ID} $i ${CLASS}_${FAMILY_ID}
    echo "bash /cluster/project1/ProCovar/colabfold-customtemplates/main_db.sh /home/dbuchan/inputs/alphafold/${CLASS}/${FASTA_ID} $i ${CLASS}_${FAMILY_ID}_${i:57:-4}"
done

# move the model
cp ${CLASS}_${FAMILY_ID}_*_unrelaxed_rank_1_model_*.pdb /home/dbuchan/alphafold_runs/
# tidy up
rm -rf ${CLASS}_${FAMILY_ID}_*_all
rm -rf ${CLASS}_${FAMILY_ID}_*_template
rm -f ${CLASS}_${FAMILY_ID}*