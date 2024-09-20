#!/bin/bash -l
#$ -o /home/ucbcdwb/Scratch/output/alpha/std.out
#$ -e /home/ucbcdwb/Scratch/output/alpha/std.err

# Request an hour of run time
#$ -l h_rt=01:00:0

# Request 10 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=10G

#$ -N runalpha
#$ -l tmem=2G
#$ -l gpu=true
#$ -pe gpu 1
#$ -R y

# Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
#$ -t 1-10000
# Set the name of the job.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/cuda-11.8/lib64/
conda activate /cluster/project1/ProCovar/lmoffat/miniconda/envs/alphafold
cd /cluster/project1/ProCovar/colabfold-customtemplates
bash /cluster/project1/ProCovar/colabfold-customtemplates/main_db.sh test_run.single_sequence.fa test_run.single_sequence.a3m db_test

# move the model
cp db_test_*_unrelaxed_rank_1_model_*.pdb /home/dbuchan/alphafold_runs/
# tidy up
rm -rf db_test_*_all
rm -rf db_test_*_template
rm -f db_test*