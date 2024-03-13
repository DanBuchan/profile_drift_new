# Drift paper

1. What does drift look like for Pfam families?
2. What does sequence generation look like for Pfam/HMMER HMMS
3. What does sequence generation look like for ESM on a single-seq and MSA basis
4. Is drift of representation-ness correlated to performance in predictive tasks, namely CASP 

## 1. Calculating PFAM distance matrices

PFAM RELEASE 36 @ 8 March 2024

1. Get the Pfam-A.full.uniprot sequences and Pfam-A.hmm model sdataset from Interpro
2. extract all pfam sequences to fasta and annotation with family ID: ~/bin/prep_pfam_fasta.py
3. Convert to PSI-BLAST DB
4. Using Hmmemit output the consensus sequence for each family: 
   ~/Applications/hmmer-3.3.2/src/hmmemit -c ~/Data/pfam/Pfam-A.hmm > pfam_consensus_reps.fa
   relabel_pfam_consensus_seqs.py ~/Data/pfam/Pfam-A.hmm ~/Data/pfam/pfam_consensus_reps.fa > pfam_consensus_reps_labelled.fa
5. Using RaxML build the distance matrix over all Pfam domains, using mafft and raxml
   calculate_pfam_distances.py

   This gives us an all-against-all evolutionary distance matrix but built on an MAFFT MSA that may not be meaningful

6. Perform an all-against-all Needleman and wunsch of the reps. Extract bits scores for a similarity matrix Scale/Normalise to between 0 and one and invert for a distance matrix
   use pfam_reps_nw.py over our relabelled pfam_consensus_reps_labelled.fa
   wrote a morcambe script run_pfam.sh to batch job this over a couple of days as it is A LOT of comparisons 20k x 20k

7. Script that combines the distances down to a big matrix, maybe numpy and save as blob or pickle.

## 1b Cluster analysis

a. Take the distance matrices and cluster them.
b. project with t-sne
c. Are the clusters meaningful?

## 2. Drift analysis

1. Psiblast Consensus seq against pfam blast db, Save number of hits for each family at each iteration, Additionally save sequences at each iteration IF we detect drift and build an MSA

## 2b 

1. How many and What families drift?
2. What kinds of drift pattern are there?
3. Do drift families share the same folds or not? If not what happens if we build alphafold models with the MSAs at each iteration (see above)
4. Do drifting families correlate with "hard" targets for alphafold/dmplfold. i.e. regions of low plDDT or just the poor models in CASP15

## 3. Other Model coherence.

1. Generate 100 seqs for each drift family using hhemit
2. blast each seq against pfam seq db. What is the family of top hit? What is the consensus family of the top 5 and 10 hits?
3. Repeat for protein ESM single and ESM-MSA


