#$ -S /bin/bash
#$ -l tmem=1G
#$ -l h_vmem=2G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N pfam_nw
#$ -t 1-496

# Run the application.
echo 'python /home/dbuchan/profile_drift_new/scripts/pfam_reps_nw.py /home/dbuchan/pfam_nw/pfam_consensus_reps_labelled_flattened.fa /home/dbuchan/pfam_nw/${SGE_TASK_ID}_pfam_consensus > /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.csv 2> /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.err'
mkdir /scratch0/pfam_nw_${SGE_TASK_ID}/
cd /scratch0/pfam_nw_${SGE_TASK_ID}/
python /home/dbuchan/profile_drift_new/scripts/pfam_reps_nw.py /home/dbuchan/pfam_nw/pfam_consensus_reps_labelled_flattened.fa /home/dbuchan/pfam_nw/${SGE_TASK_ID}_pfam_consensus > /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.csv 2> /home/dbuchan/pfam_nw/${SGE_TASK_ID}_hits.err
cd ~/
rm -rf /scratch0/pfam_nw_${SGE_TASK_ID}
