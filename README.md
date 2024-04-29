# Drift paper

1. What does drift look like for Pfam families?
2. What does sequence generation look like for Pfam/HMMER HMMS
3. What does sequence generation look like for ESM on a single-seq and MSA basis
4. Is drift of representation-ness correlated to performance in predictive tasks, namely CASP 

## 1. Calculating PFAM distance matrices

PFAM RELEASE 36 @ 15 March 2024

1. Get the Pfam-A.full.uniprot sequences and Pfam-A.hmm model sdataset from Interpro
2. extract all pfam sequences to fasta and annotation with family ID: ~/bin/prep_pfam_fasta.py
   discarding truncated sequence
3. Convert to PSI-BLAST DB
4. Using Hmmemit output the consensus sequence for each family: 
   > ~/Applications/hmmer-3.3.2/src/hmmemit -c ~/Data/pfam/Pfam-A.hmm > pfam_consensus_reps.fa
   > python relabel_pfam_consensus_seqs.py ~/Data/pfam/Pfam-A.hmm ~/Data/pfam/pfam_consensus_reps.fa > pfam_consensus_reps_labelled.fa
   > python ~/bin/prep_fasta.py pfam_consensus_reps_labelled.fa > pfam_consensus_reps_labelled_flattened.fa
5. Using RaxML build the distance matrix over all Pfam domains, using mafft and raxml
   >python3 ./calculate_pfam_distances.py ~/Data/pfam/pfam_consensus_reps_labelled_flattened.fa
   This gives us an all-against-all evolutionary distance matrix but built on an MAFFT MSA that may not be meaningful. But fun to do/try

6. Perform an all-against-all Needleman and wunsch of the reps. Extract bits scores for a similarity matrix Scale/Normalise to between 0 and one and invert for a distance matrix
   > wc -l pfam_consensus_reps_labelled_flattened.fa
   Then use split to divide in to 500 files (880 will change depending on size), for cluster execution
   > split --numeric-suffixes=1 -a 3 -l 84 --additional-suffix=_pfam_consensus pfam_consensus_reps_labelled_flattened.fa ''

   use pfam_reps_nw.py over our relabelled pfam_consensus_reps_labelled.fa
   wrote a morecambe script run_pfam_nw.sh to batch job this over a couple of days as it is A LOT of comparisons 20k x 20k

7. Script that combines the distances down to a big matrix, maybe numpy and save as blob or pickle.
   build_pfam_distance_matrix.py

## 1b Cluster analysis

a. Take the distance matrices and cluster them.
b. project with t-sne
c. Are the clusters meaningful?

## 2. Drift analysis

1. Use prep_pfam_fasta.py to make a fasta file of all PF families in Pfam-A.full.uniprot, The makeblastdb to make a blsat db of it
2. Psiblast Consensus seq against pfam blast db, Save number of hits for each family at each iteration, Additionally save sequences at each iteration IF we detect drift and build an MSA
   > python run_pfam_rep_blasts.py
  use run_pfam_psiblast.sh on myriad to coordinate this.
3. Analyse *_blast_summary.csv, to find out which families show drift and what kinds of drift
   > python calculate_drift_types.py results_data/psiblast_iteration_summaries/
   outputs
   drift_list.txt - familes with drift contamination
   non_drift_list.txt - families with no such contamination
   drift_error_list.txt - families with errors in the drift summary file, to be re-run
   insignificant_drifts.txt - families where contamination is less than 5% of the peak family size
   significant_drifts.txt - list of drift families with substantial amounts of drift
   drift_summary.txt - summary of some different behaviours of the query and contaminant families
4. We select some Pfam families from each drift class. Unpack the sequences. Align them with mafft.
   > select_targets.py
   we take in some lists of families and output a list of targets that we want to model
   alpha_fold_targets.csv - a list of target pfam families for differing drift classes
   > get_target_id_list.py
   little helper script that takes the target list and translates it in to the JOB IDs that were used on the cluster to output the target list
   mafft_targets.txt - just the list of JOB IDs that map to the alphafold targets
   > run_mafft.py mafft.sh
5. Take all our alignments and then run alphafold2.
6. Collate models and analyse.


# 3. Model coherence

1. Generate 100 seqs for each hhemit, take in alpha_fold_targets/mafft_targets.txt and output 100 seqs for each drift class
   > hmmer_seq_generator.py
   hmm_subset.hmm - the hmms of the target subset
   hmm_generated_seqs.fa - 100 seqs per family in the target list
2. Find the pfam family that our mafft/alpha fold targets best hit using fasta
   find_closest_hmm_seqfamily.py run_closest.hmm.sh
   outputs lots of ".best" files that we cat in to one large csv


## 2b 

1. Analyse *_blast_summary.csv, to find out which families show drift and what kinds of drift

1. How many and What families drift?
2. What kinds of drift pattern are there?
3. Do drift families share the same folds or not? If not what happens if we build alphafold models with the MSAs at each iteration (see above)
4. Do drifting families correlate with "hard" targets for alphafold/dmplfold. i.e. regions of low plDDT or just the poor models in CASP15

## 3. Other Model coherence.

1. Generate 100 seqs for each drift family using hhemit
2. blast each seq against pfam seq db. What is the family of top hit? What is the consensus family of the top 5 and 10 hits?
3. Repeat for protein ESM single and ESM-MSA


