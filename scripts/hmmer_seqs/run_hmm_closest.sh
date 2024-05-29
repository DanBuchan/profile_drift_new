#$ -S /bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=4:0:0
#$ -j y
#$ -N hmm_closest
##$ -t 1-51900
#$ -t 10001-20000

# Run the application.

mkdir /scratch0/pfam_nw_${SGE_TASK_ID}/
cd /scratch0/pfam_nw_${SGE_TASK_ID}/

echo 'python scripts/hmmer_seqs/find_closest_hmm_seq_family.py /home/dbuchan/inputs/hmm_generated_seqs.fa /home/dbuchan/inputs/Pfam-A.full.uniprot.fa /home/dbuchan/outputs/hmm_closest/ ${SGE_TASK_ID}'
python /home/dbuchan/profile_drift_new/scripts/hmmer_seqs/find_closest_hmm_seq_family.py /home/dbuchan/inputs/hmm_generated_seqs.fa /home/dbuchan/inputs/Pfam-A.full.uniprot.fa /home/dbuchan/outputs/hmm_closest/ ${SGE_TASK_ID}

cd ~/
rm -rf /scratch0/pfam_nw_${SGE_TASK_ID}
