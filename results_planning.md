1. Chart/table of how many drift and proportions of the types of drift
    - pfam_rep_drift_summary/
    20,049 pfam families
    14,273 - show no signs of drift (71%)
    5,776 - find sequence from more than one pfam family
           - of these 2,748 never have more than 5% of the sequences as contaminants (48%)
           - 3,028 (52%) have significant contamination from non-family domains (15% of all the searches)

    then a table of drift_overview.txt

From that we sample 530 families

2. Chart table of hmm seqs, how many are in family and how many are out of family
    - /results_data/drift/best_hits/hmm_closest/drift_summary.csv 
    - results_data/generation_or_af_targets/alphafold_targets.csv
    - hmm_drift_percentages.csv
3. Chart table of prottrans, how many are in family and how many are out of family, broken down by 25, 50, 75 mask
- /results_data/drift/best_hits/prottrans_closest/drift_summary.csv
- prottrans_drift_percentages.csv


4. Model analysis
  - results_data/generation_or_af_targets/alphafold_targets.csv