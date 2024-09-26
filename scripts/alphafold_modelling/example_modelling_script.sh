#!/bin/bash -l

#$ -S /bin/bash
# Request an hour of run time
#$ -l h_rt=300:00:0

# Request 10 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tscratch=1G
#$ -l scratch0free=1G

#$ -N runalpha
#$ -l tmem=2G
#$ -l gpu=true
#$ -pe gpu 1
#$ -R y

# Set up the job array.  In this instance we have requested 10000 tasks
# numbered 1 to 10000.
# Set the name of the job.

#$ -t 1-6

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/cuda-11.8/lib64/
/home/dbuchan/miniconda3/condabin/conda init bash
/home/dbuchan/miniconda3/condabin/conda activate /cluster/project1/ProCovar/lmoffat/miniconda/envs/alphafold 
cd /cluster/project1/ProCovar/colabfold-customtemplates

CLASSES=('contaminants_grew' 'contaminants_complex' 'contaminants_purified' 'insig_drift' 'non_drift' 'query_purified')
CLASS=${CLASSES[${SGE_TASK_ID}]}

FAMILY_ID=${FASTA_ID::-3}

for i in /home/dbuchan/inputs/alphafold/${CLASS}/*.a2m; do
    FAMILY_ID=`echo ${i} | grep -oP "PF\d{5}"`
    ENDING_CONFIG=`echo ${i} | grep -oP "PF\d{5}_\d+"`
    # bash /cluster/project1/ProCovar/colabfold-customtemplates/main_db.sh /home/dbuchan/inputs/alphafold/${CLASS}/${FASTA_ID}.fa $i ${CLASS}_${ENDING_CONFIG}
    echo "bash /cluster/project1/ProCovar/colabfold-customtemplates/main_db.sh /home/dbuchan/inputs/alphafold/${CLASS}/${FASTA_ID} $i ${CLASS}_${ENDING_CONFIG}"
done

# move the model
cp ${CLASS}_${FAMILY_ID}_*_unrelaxed_rank_1_model_*.pdb /home/dbuchan/alphafold_runs/
# tidy up
rm -rf ${CLASS}_${FAMILY_ID}_*_all
rm -rf ${CLASS}_${FAMILY_ID}_*_template
rm -f ${CLASS}_${FAMILY_ID}*