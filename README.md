# Drift paper

1. What does drift look like for Pfam families?
2. What does sequence generation look like for Pfam/HMMER HMMS
3. What does sequence generation look like for ESM on a single-seq and MSA basis
4. Is drift of representation-ness correlated to performance in predictive tasks, namely CASP 

## 1. Calculating drift in Pfam families

PFAM RELEASE 36 @ 8 March 2024

a. Get the Pfam-A.full.uniprot sequences and Pfam-A.hmm model sdataset from Interpro
b. extract all pfam sequences to fasta and annotation with family ID: ~/bin/prep_pfam_fasta.py
c. Convert to PSI-BLAST DB
d. Using Hmmemit output the consensus sequence for each family: 
   ~/Applications/hmmer-3.3.2/src/hmmemit -c ~/Data/pfam/Pfam-A.hmm > pfam_consensus_reps.fa
   relabel_pfam_consensus_seqs.py ~/Data/pfam/Pfam-A.hmm ~/Data/pfam/pfam_consensus_reps.fa > pfam_consensus_reps_labelled.fa
e. Using RaxML build the distance matrix over all Pfam domains, using mafft and raxml
   calculate_pfam_distances.py

   This gives us an all-against-all evoltionary distance matrix but built on an MAFFT MSA that may not be meaningful

f. Perform an all-against-all Needleman and wunsch of the reps. Extract bits scores for a similarity matrix Scale/Normalise to between 0 and one and invert for a distance matrix
   use pfam_reps_nw.py over our relabelled pfam_consensus_reps_labelled.fa
   wrote a morcambe script run_pfam.sh to batch job this over a couple of days as it is A LOT of comparisons 20k x 20k


g. Psiblast Consensus seq against pfam blast db
h. record drift

## Cluster analysis

a. Take the distance matrices and cluster them.
b. project with t-sne
c. Are the clusters meaningful?

## 2. For 
