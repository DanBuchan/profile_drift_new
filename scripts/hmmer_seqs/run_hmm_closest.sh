#!/bin/bash -l
#$ -o /home/ucbcdwb/Scratch/output/profile/std.out
#$ -e /home/ucbcdwb/Scratch/output/profile/std.err

# Request an hour of run time
#$ -l h_rt=01:00:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=8G

# Request 10 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=8G

# Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
#$ -t 1-26500

# Set the name of the job.

#$ -N hmmGene

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/ucbcdwb/Scratch/output/hmm/

mkdir /home/ucbcdwb/Scratch/output/hmm/hmm_fa_${SGE_TASK_ID}/
cd /home/ucbcdwb/Scratch/output/hmm/hmm_fa_${SGE_TASK_ID}/

echo 'python scripts/hmmer_seqs/find_closest_hmm_seq_family.py /home/dbuchan/inputs/hmm_generated_seqs.fa /home/dbuchan/inputs/Pfam-A.full.uniprot.fa /home/dbuchan/outputs/hmm_closest/ ${SGE_TASK_ID}'
python /home/ucbcdwb/Applications/profile_drift_new/scripts/hmmer_seqs/find_closest_hmm_seq_family.py /home/ucbcdwb/inputs/hmm_generated_seqs.fa /home/ucbcdwb/inputs/Pfam-A.full.uniprot.fa /home/ucbcdwb/Scratch/output/hmm ${SGE_TASK_ID}

cd ~/
rm -rf /home/ucbcdwb/Scratch/output/hmm/hmm_fa_${SGE_TASK_ID}
