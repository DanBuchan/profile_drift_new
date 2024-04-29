#!/bin/bash -l
#$ -o /home/ucbcdwb/Scratch/output/profile/std.out
#$ -e /home/ucbcdwb/Scratch/output/profile/std.err

# Request an hour of run time
#$ -l h_rt=60:00:0

#$ -l mem=16G

# Request 10 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=10G

# Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
#$ -t 1-519

# Set the name of the job.

#$ -N runMAFFT

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/ucbcdwb/Scratch/output/profile/

# Run the application.
module load python3
source /home/ucbcdwb/Scratch/virtualenvs/profile_drift/bin/activate
cd /home/ucbcdwb/Scratch/output/profile/
echo "python /home/ucbcdwb/Applications/profile_drift_new/scripts/run_mafft.py $SGE_TASK_ID /home/ucbcdwb/Applications/profile_drift_new//mafft_targets.txt"
python /home/ucbcdwb/Applications/profile_drift_new/scripts/run_mafft.py $SGE_TASK_ID /home/ucbcdwb/Applications/profile_drift_new//mafft_targets.txt
