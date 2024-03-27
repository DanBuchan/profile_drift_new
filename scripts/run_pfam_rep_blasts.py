from collections import defaultdict
import subprocess
import sys
import os
import re
import glob
from Bio.Blast import NCBIXML

# python run_pfam_rep_blasts.py ~/Data/pfam/test.fa ~/Data/pfam/Pfam-A.full.uniprot.fa [ID NUMBER [ITERATIONS]

def read_reps(dom_seqs):
    """
    read the domain seqs and pull out a random rep for each. Not really caring
    about the cath assigned rep as that's somewhat arbitrary anyway
    """
    reps = defaultdict(dict)
    rep_id = ''
    family = ''
    with open(dom_seqs, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                line = line.rstrip()
                entries = line.split("|")
                family = entries[1]
                rep_id = entries[0].split("/")[0][1:]
            else:
                seq = line.rstrip()
                reps[family] = {rep_id: seq}
    return reps

def do_blast_iterations(id, blast_db, n):
    in_file = f'{id}.fa'
    #first_iteration_args = ['/home/dbuchan/Applications/ncbi-blast-2.12.0+/bin/psiblast',
    first_iteration_args = ['/home/ucbcdwb/Applications/ncbi-blast-2.12.0+/bin/psiblast',
                            '-query',
                            in_file,
                            '-num_iterations',
                            '1',
                            '-db',
                            blast_db,
                            '-out_pssm',
                            f'{id}_iteration1.pssm',
                            '-outfmt',
                            '5',
                            '-out',
                            f'{id}_iteration1.xml',
                            '-save_pssm_after_last_round',
                            '-max_target_seqs',
                            '50000']
    print(" ".join(first_iteration_args))
    subprocess.call(first_iteration_args)
    for i in range(2, int(n)+1):
    #    iteration_args = ['/home/dbuchan/Applications/ncbi-blast-2.12.0+/bin/psiblast',
        iteration_args = ['/home/ucbcdwb/Applications/ncbi-blast-2.12.0+/bin/psiblast',
                                '-in_pssm',
                                f'{id}_iteration{i-1}.pssm',
                                '-num_iterations',
                                '1',
                                '-db',
                                blast_db,
                                '-out_pssm',
                                f'{id}_iteration{i}.pssm',
                                '-outfmt',
                                '5',
                                '-out',
                                f'{id}_iteration{i}.xml',
                                '-save_pssm_after_last_round',
                                '-max_target_seqs',
                                '50000']
        print(" ".join(iteration_args))    
        subprocess.call(iteration_args)

def process_blast_results(file_id, seq, family, iterations):

    fhsummary = open(f"{file_id}_blast_summary.csv", "w")
    fhsummary.write("iteration,query,query_family,hit_family,count\n")
    pf_pattern = re.compile(r"\|(PF\d{5})")
    iteration_counts = {}
    for i in range(1, iterations+1):
        print("parsing", i)
        counts = defaultdict(int)
        blastfh = open(f"{file_id}_iteration{i}.xml", "r")
        fhseqsout = open(f"{file_id}_iteration{i}_seqs.fa", "w")
        fhseqsout.write(f">{file_id}\n")
        fhseqsout.write(f"{seq}\n")
        results = NCBIXML.parse(blastfh)
        for hit in results:
            for alignment in hit.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.00001:
                        fhseqsout.write(f">{alignment.title}\n")
                        fhseqsout.write(f"{hsp.sbjct}\n")

                        result = re.search(pf_pattern, alignment.title)
                        if result:
                            # print(result.groups()[0])
                            counts[result.groups()[0]] += 1
        iteration_counts[i] = counts
        blastfh.close()
        fhseqsout.close()
    # print(iteration_counts)
    for iteration in iteration_counts.keys():
        for qfamily in iteration_counts[iteration].keys():
            fhsummary.write(f"{iteration},{family},{qfamily},{iteration_counts[iteration][qfamily]}\n")
    fhsummary.close()
        

def align_seqs(seq_file, iterations):
    clean_up = False
    for i in range(1, iterations+1):
        seqs = f"{seq_file}_iteration{i}_seqs.fa"
        msa = f"{seq_file}_iteration{i}_seqs.msa"
        mafft_args = [
        #    '/home/dbuchan/Applications/mafft-7.490-with-extensions/core/mafft',
            '/home/ucbcdwb/Applications/mafft-linux64/mafft.bat',
            seqs,
        ]
        # print(" ".join(mafft_args))
        try:
            msa_data = subprocess.check_output(mafft_args)
            fhOut = open(msa, "wb")
            fhOut.write(msa_data)
            fhOut.close()
            clean_up=True
        except Exception as e:
            clean_up = False
            break
    
    if clean_up:
        tar_args = [
            '/usr/bin/tar',
            'czf',
            f'{seq_file}_msa.tar.gz',
        ]
        tar_args.extend(glob.glob("*.msa"))
        tar_output = subprocess.check_output(tar_args)
        for seqs in glob.glob("*_seqs.fa"):
            os.remove(seqs)
        for msa in glob.glob("*_seqs.msa"):
            os.remove(msa)

def run_blasts(family, id, seq, blast_db, iterations):
    """
    run psiblast over each rep against the cath dom seqs db
    """
    fhRep = open(f'{id}.fa', "w")
    fhRep.write(f">{id}|{family}\n")
    fhRep.write(f"{seq}\n")
    fhRep.close()
    # do_blast_iterations(id, blast_db, iterations)
    # process_blast_results(id, seq, family, iterations)
    align_seqs(id, iterations)
    os.remove(f'{id}.fa')
    for xml in glob.glob("*.xml"):
        os.remove(xml)
    for pssm in glob.glob("*.pssm"):
        os.remove(pssm)

def read_fasta(file):
    seqs = {}
    fasta_name = None
    seq = None
    with open(file, "r") as fhIn:
        for line in fhIn:
            if line.startswith(">"):
                if fasta_name:
                    seqs[fasta_name] = seq
                fasta_name = line.rstrip()
                seq = ''
            else:
                seq = seq+line.rstrip()
    return(seqs)

rep_file = sys.argv[1]
blast_db = sys.argv[2]

rep_seqs = read_fasta(rep_file)
id_list = list(rep_seqs.keys())

id = id_list[int(sys.argv[3])-1]
family = id[-7:]
seq_id = id[1:-8]
seq_id = seq_id.replace("/", "_")
run_blasts(family, seq_id, rep_seqs[id], blast_db, int(sys.argv[4]))

# for id in rep_seqs.keys():
#     family = id[-7:]
#     seq_id = id[1:-8]
#     seq_id = seq_id.replace("/", "_")
#     run_blasts(family, seq_id, rep_seqs[id], blast_db)
#     exit()