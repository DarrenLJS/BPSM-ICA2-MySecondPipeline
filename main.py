#!.venv/bin/python3

import json
import subprocess
import re
import sys
from input_handler import validate_protein, validate_taxon
from fetch_sequence import run_esearch, run_efetch
from parse_fasta import parse_fasta
from conservation_analysis import temp_fasta, run_clustalo, run_plotcon
from scan_prosite import run_prosite_scan, parse_prosite_output
from blast_analysis import make_blast_db, select_blast_ref, run_blast_search, parse_blast_output

def main():

    while True:
        try:
            pfam_input = input("Protein family: ").strip()
            pfam_id, pfam_name = validate_protein(pfam_input)
            break
        except Exception as e:
            print(f"\nPlease input a valid protein family. {e}")

    while True:
        email = input("Your email (needed for NCBI URL request): ").strip()
        if re.search(
            r"^[a-zA-Z0-9.!#$%&'*+\/=?^_`{|}~-]+@[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?(?:\.[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?)*$", 
            email, re.IGNORECASE
        ):
            break
        else:
            print("\nPlease input a valid email.")

    while True:
        try:
            taxon_input = input("Taxon group: ").strip()
            taxon_id, taxon_name = validate_taxon(taxon_input, email)
            break
        except Exception as e:
            print(f"\nPlease input a valid taxon group. {e}")
    
    print(f"\nSuccess! Your protein family: {pfam_name}\tPfam ID: {pfam_id}")
    print(f"Success! Your taxon group: {taxon_name}\tTaxon ID: {taxon_id}\n")
    
    print("Searching NCBI Protein...")
    ids_list = run_esearch(pfam_name, taxon_id, email)
    if not ids_list:
        sys.exit("\nNo UIDs could be retrieved from NCBI Protein! Please try again.")

    elif len(ids_list) > 1000:
        while True:
            try:
                warner = input("Your dataset contains more than 1000 protein sequences!\nAre you sure you want to continue? (y/n) ").strip().lower()
                if warner not in ["y", "n", "yes", "no"]:
                    raise ValueError
                break
            except Exception as e:
                print("\nInvalid input.")
        if warner in ["n", "no"]:
            sys.exit("\nPipeline stopped!")
    
    print("Protein UIDs successfully retrieved!")
    out_dir = f"{pfam_name}_{taxon_name}"
    subprocess.call(f"rm -rf {out_dir}", shell = True)
    subprocess.call(f"mkdir -p {out_dir}", shell = True)
    
    print("Fetching all protein sequences in FASTA format...")
    out_fasta = run_efetch(ids_list, out_dir, email)
    print(f"FASTA sequences successfully downloaded to {out_dir}/{out_fasta}")

    print(f"Parsing {out_dir}/{out_fasta}...")
    records = parse_fasta(out_dir, out_fasta)
    with open("test_json.json", "w") as file:
        json.dump(records, file, indent = 2)
    print(f"Successfully parsed {out_dir}/{out_fasta}")
    
    #out_temp_fasta = temp_fasta(records, out_dir)
    #run_clustalo(out_dir, out_temp_fasta)
    #run_plotcon(out_dir)

    #out_motifs = run_prosite_scan(out_dir, out_fetch)
    #parse_prosite_output(out_dir, out_motifs)
    
    #blastdb_name = f"{out_dir}_blastdb"
    #make_blast_db(out_dir, out_fetch, blastdb_name)
    #out_reference = select_blast_ref(records, out_dir)
    #out_blast = run_blast_search(out_dir, out_reference, blastdb_name)
    #parse_blast_output(out_dir, out_blast)


if __name__ == "__main__":
    main()
