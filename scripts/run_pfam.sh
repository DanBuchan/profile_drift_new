#$ -S /bin/bash
#$ -l tmem=1G
#$ -l h_vmem=2G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N pfam_nw
#$ -t 1-50

# Run the application.
/home/dbuchan/profile_drift_new/scripts/pfam_reps_nw.py /home/dbuchan/pfam_nw/pfam_consensus_reps_labelled_flattened.fa /home/dbuchan/pfam_nw/${SGE_TASK_ID}_pfam_consensus > /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.csv 2> /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.err