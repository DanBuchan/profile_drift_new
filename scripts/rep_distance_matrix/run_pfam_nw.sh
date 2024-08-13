#$ -S /bin/bash
#$ -l tmem=1G
#$ -l tscratch=1G
#$ -l h_vmem=2G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N pfam_nw
#$ -t 3-504

# Run the application.
hostname
mkdir /scratch0/pfam_nw_${SGE_TASK_ID}/
cd /scratch0/pfam_nw_${SGE_TASK_ID}/
cp /home/dbuchan/pfam_nw/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/
cp /home/dbuchan/pfam_nw/${SGE_TASK_ID}_pfam_random /scratch0/pfam_nw_${SGE_TASK_ID}/
echo 'python /home/dbuchan/profile_drift_new/scripts/rep_distance_matrix/pfam_reps_nw.py /scratch0/pfam_nw_${SGE_TASK_ID}/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_pfam_random > /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.csv 2> /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.err'
python /home/dbuchan/profile_drift_new/scripts/rep_distance_matrix/pfam_reps_nw.py /scratch0/pfam_nw_${SGE_TASK_ID}/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_pfam_random > /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.csv 2> /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.err
cp /scratch0/pfam_nw_${SGE_TASK_ID}/*.csv /home/dbuchan/pfam_nw/
cp /scratch0/pfam_nw_${SGE_TASK_ID}/*.err /home/dbuchan/pfam_nw/
cd ~/
rm -rf /scratch0/pfam_nw_${SGE_TASK_ID}
