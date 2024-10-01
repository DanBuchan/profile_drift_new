# Models to build

1. non_drift - 2 models, 1st and the last : 100 examples == 200 models
2. insig_drift - 2 models build the 1st and the last : 100 examples == 200 models
3. contaminants_grew - 4 models, 1st, last and 5th and 15th : 100 examples == 400 models
4. contaminants_purified - 3 models 1st, last and peak contamination : 50 examples = 150 models
5. contaminants_complex -  3 models 1st, last and peak contamination : 100 examples = 300 models
6. query_purified - 4 models, 1st, last and 5th and 15th : 100 examples == 400 models

total 1650 models == 23 GPU days

Or we could just do all 11,000 models == 153 GPU days