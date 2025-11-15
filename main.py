#!/usr/bin/python3

import subprocess
from input_handler import get_input, validate_input
from fetch_sequence import run_esearch, run_efetch
from parse_fasta import parse_fasta
from conservation_analysis import temp_fasta, run_clustalo, run_plotcon
from scan_prosite import run_prosite_scan, parse_prosite_output
from blast_analysis import make_blast_db, select_blast_ref, run_blast_search, parse_blast_output

def main():

    while True:
        try:
            protein_family, taxon_group = get_input()
            validate_input(protein_family, taxon_group)
            break
        except ValueError:
            pass
    
    out_dir = f"{protein_family}_{taxon_group}"
    subprocess.call(f"rm -rf {out_dir}")
    subprocess.call(f"mkdir -p {out_dir}")
    
    ids_list = run_esearch(protein_family, taxon_group)
    out_fetch = run_efetch(ids_list, out_dir)
    
    records = parse_fasta(out_dir, out_fetch)
    
    out_temp_fasta = temp_fasta(records, out_dir)
    run_clustalo(out_dir, out_temp_fasta)
    run_plotcon(out_dir)

    out_motifs = run_prosite_scan(out_dir, out_fetch)
    parse_prosite_output(out_dir, out_motifs)
    
    blastdb_name = f"{out_dir}_blastdb"
    make_blast_db(out_dir, out_fetch, blastdb_name)
    out_reference = select_blast_ref(records, out_dir)
    out_blast = run_blast_search(out_dir, out_reference, blastdb_name)
    parse_blast_output(out_dir, out_blast)


if __name__ == "__main__":
    main()
