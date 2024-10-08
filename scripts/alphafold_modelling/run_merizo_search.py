import glob
import sys
import re
from subprocess import Popen, PIPE

# Run merizo search cath over all our alphafold models to see what H family they are
# python run_merizo_search.py ./results_data/alphafold_models/
#


# python ./merizo_search/merizo.py search ~/Projects/profile_drift/results_data/alphafold_models/query_purified_PF00050_1_unrelaxed_rank_1_model.pdb ./examples/database/cath query_purified_PF00050_1 tmp --format query,emb_rank,target,emb_score,q_len,t_len,ali_len,seq_id,q_tm,t_tm,max_tm,rmsd,metadata --output_headers


for file in glob.glob(f'{sys.argv[1]}/*.pdb'):
    m = re.search('models/(.+_PF\d+_\d+)_', file)
    identifier = ''
    if m:
        identifier = m.group(1)
    # print(file)
    if "contaminants_complex_PF00106_20_unrelaxed_rank_1_model.pdb" not in file:
        continue
    # print(m.group(1))
    print(identifier)
    args = ['python',
            '/home/dbuchan/Code/merizo_search/merizo_search/merizo.py',
            'search',
            file,
            '/home/dbuchan/Code/merizo_search/examples/database/cath',
            identifier,
            'tmp',
            '--format',
            'query,emb_rank,target,emb_score,q_len,t_len,ali_len,seq_id,q_tm,t_tm,max_tm,rmsd,metadata',
            '--output_headers',
    ]
    print("Calculating", " ".join(args))
    try:
        p = Popen(args, stdout=PIPE, stderr=PIPE)
        result_stdout, err = p.communicate()
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if p.returncode != 0:
        print("Non Zero Exit status: "+str(p.returncode))
        raise OSError("Non Zero Exit status: "+str(p.returncode))
    results_file = f'{identifier}_search.tsv'
    print(results_file)
    with open(results_file, "r", encoding="utf-8") as fhIn:
        for line in fhIn:
            entries = line.split("\t")
            hit = entries[2]
            max_tm = entries[10]
            meta = entries[12][1:-2]
            meta_fields = meta.split(", ")
            print(meta_fields)
            cath_fields = meta_fields[0].split(": ")
            print(cath_fields)
            h_family = cath_fields[1]
            print(hit, max_tm, cath_fields)
            
    # break