#!.venv/bin/python3

import subprocess

def make_blast_db(out_dir, fasta_filename, blast_db):
    subprocess.call(f"makeblastdb -in {out_dir}/{fasta_filename} -dbtype prot -out {out_dir}/{out_dir}_blast", shell = True)

def select_blast_ref(records, out_dir):
    return out_filename

def run_blast_search(out_dir, blastref_filename, blast_db):
    return out_filename

def parse_blast_output(out_dir, blastresults_filename):
    pass

