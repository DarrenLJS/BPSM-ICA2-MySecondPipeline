#!.venv/bin/python3

import subprocess

def make_blast_db(out_dir, fasta_filename):
    blast_db = f"{out_dir}_blast"
    subprocess.call(f"makeblastdb -in {out_dir}/{fasta_filename} -dbtype prot -out {out_dir}/{blast_db}", shell = True)
    return blast_db

def select_blast_ref(records, out_dir):
    ref_seq = list(sorted(records, key = lambda x: x["length"], reverse = True))[0]
    blast_ref = f"{out_dir}_blast_ref.fasta"
    with open(f"{out_dir}/{blast_ref}", "w") as file:
        file.write(f">{ref_seq['uid']} {ref_seq['description']} [{ref_seq['species']}]\n")
        file.write(f"{ref_seq['sequence']}\n\n")
    return blast_ref

def run_blast_search(out_dir, blastref_filename, blast_db):
    blast_output = f"{out_dir}_blastoutput.out"
    subprocess.call(f"blastp -db {out_dir}/{blast_db} -query {out_dir}/{blastref_filename} -outfmt 7 > {out_dir}/{blast_output}", shell = True)
    return blast_output

def parse_blast_output(out_dir, blastresults_filename):
    pass

