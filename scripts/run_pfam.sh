#$ -S /bin/bash
#$ -l tmem=1G
#$ -l h_vmem=2G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N pfam_nw

# Run the application.
module load python3
source /home/ucbcdwb/Scratch/profile_drift/profile_drift/bin/activate
mkdir $SGE_TASK_ID
cd $SGE_TASK_ID/
echo "python run_pfam_rep_blasts.py ~/Scratch/Data/pfam/reps.fasta.fa ~/Scratch/Data/pfam/pfam_fasta.fa $SGE_TASK_ID"
python /home/ucbcdwb/Scratch/profile_drift/run_pfam_rep_blasts.py ~/Scratch/Data/pfam/reps.fasta.fa ~/Scratch/Data/pfam/pfam_fasta.fa $SGE_TASK_ID
