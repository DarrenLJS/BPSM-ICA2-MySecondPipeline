#!.venv/bin/python3

import subprocess

def build_prosite_input(records, out_dir):
    prosite_input = f"{out_dir}_prosite_in.fasta"
    with open(f"{out_dir}/{prosite_input}", "w") as file:
        for entry in records:
            file.write(f">{entry['uid']} {entry['description']} [{entry['species']}]\n")
            file.write(f"{entry['sequence']}\n\n")
    return prosite_input

def run_prosite_scan(out_dir, prosite_input_filename):
    prosite_output = f"{out_dir}_prosite_out.patmatmotifs"
    subprocess.call(f"patmatmotifs -sequence {out_dir}/{prosite_input_filename} -outfile {out_dir}/{prosite_output}", shell = True)
    return prosite_output

def parse_prosite_output(out_dir, motifs_filename):
    pass

