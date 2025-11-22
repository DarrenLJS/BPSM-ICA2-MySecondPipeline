#!.venv/bin/python3

import glob
import json
import subprocess
import re
import sys
import pandas as pd
from input_handler import validate_protein, validate_taxon
from fetch_sequence import run_esearch, run_efetch
from parse_fasta import parse_fasta
from conservation_analysis import build_clustalo_input, run_clustalo, run_plotcon, get_consensus, get_infoalign
from scan_prosite import parse_clu_results, build_prosite_input, run_prosite_scan, parse_prosite_output
from blast_analysis import make_blast_db, select_blast_ref, run_blast_search, parse_blast_output

def main():
    
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
            pfam_input = input("Protein family: ").strip()
            pfam_id, pfam_name = validate_protein(pfam_input)
            break
        except Exception as e:
            print(f"\nPlease input a valid protein family. {e}")

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
    
    print("Protein UIDs successfully retrieved!\n")
    out_dir = f"{pfam_name}_{taxon_name}"
    subprocess.call(f"rm -rf {out_dir}", shell = True)
    subprocess.call(f"mkdir -p {out_dir}", shell = True)
    
    print("Fetching all protein sequences in FASTA format...")
    out_fasta = run_efetch(ids_list, out_dir, email)
    print(f"FASTA sequences successfully downloaded to {out_fasta}\n")

    print(f"Parsing {out_fasta}...")
    records = parse_fasta(out_dir, out_fasta)
    with open(f"{out_dir}/{out_dir}_records.json", "w") as file:
        json.dump(records, file, indent = 2)
    print(f"Successfully parsed {out_fasta}\n")
    
    while True:
        try:
            print(f"Please input a number for the number of sequences to be analysed with Clustalo. Maximum: {len(records)}")
            print("The top few sequences with the largest lengths will be used for conservation analysis.")
            print("Note: The larger the number, the slower the analysis. Recommended: 1000-2000")
            conservation_analysis_size = int(input(f"Sample size: ").strip())
            if conservation_analysis_size > len(records):
                raise ValueError("Invalid input exceeds maximum size.")

            if conservation_analysis_size > 2000:
                while True:
                    try:
                        warner2 = input("\nYour input exceeds the recommended range. Are you sure you want to continue? (y/n) ").strip().lower()
                        if warner2 not in ["y", "n", "yes", "no"]:
                            raise ValueError
                        break
                    except Exception as e:
                        print("\nInvalid input.")
                if warner2 in ["n", "no"]:
                    print("\nRetrying...")
                    continue
            
            break
        except Exception as e:
            print(f"\nInvalid input. {e}")
    
    print(f"\nYour sample size: {conservation_analysis_size}")
    print(f"The {conservation_analysis_size} longest protein sequences will be used for conservation analysis.\n")

    print("Performing conservation analysis with ClustalO and Plotcon...")
    out_clustalo = f"{out_dir}_clustalo"
    subprocess.call(f"mkdir -p {out_dir}/{out_clustalo}", shell = True)
    clustalo_input = build_clustalo_input(records, out_dir, out_clustalo, conservation_analysis_size)
    clustalo_output = run_clustalo(out_dir, out_clustalo, clustalo_input)
    print(f"ClustalO output successfully generated! Generated {clustalo_output}")
    plotcon_graph, plotcon_text = run_plotcon(out_dir, out_clustalo, clustalo_output)
    print(f"Plotcon output successfully generated! Generated {plotcon_graph}, {plotcon_text}")
    consensus = get_consensus(out_dir, out_clustalo, clustalo_output)
    print(f"Successfully obtained consensus sequence! Generated {consensus}")
    infoalign_results = get_infoalign(out_dir, out_clustalo, clustalo_output)
    print(f"Successfully obtained infoalign results for {clustalo_output}! Generated {infoalign_results}\n")
    
    while True:
        try:
            print(f"Please input a number for the number of sequences to be scanned for PROSITE motifs. Maximum: {conservation_analysis_size}")
            print("The top few sequences with the highest similarity with the consensus sequence will be scanned for PROSITE motifs.")
            scan_size = int(input(f"Scan size: ").strip())
            if scan_size > conservation_analysis_size:
                raise ValueError("Invalid input exceeds maximum size.")
            break
        except Exception as e:
            print(f"\nInvalid input. {e}")
    
    print(f"\nYour scan size: {scan_size}")
    print(f"The {scan_size} protein sequences with the highest similarity to the consensus sequence will be scanned for PROSITE motifs.\n")

    out_prosite = f"{out_dir}_prosite"
    subprocess.call(f"mkdir -p {out_dir}/{out_prosite}", shell = True)
    top_records = parse_clu_results(records, out_dir, out_clustalo, infoalign_results, out_prosite, scan_size)
    print("Scanning protein sequences for PROSITE motifs...")
    prosite_input = build_prosite_input(top_records, out_dir, out_prosite)
    run_prosite_scan(out_dir, out_prosite, prosite_input)
    prosite_output = list(map(lambda x: x + ".patmatmotifs", prosite_input))
    print(f"PROSITE output successfully generated! Generated {len(prosite_output)} .patmatmotifs files.")
    prosite_summary_df, prosite_motifs = parse_prosite_output(out_dir, out_prosite, prosite_output)
    print(f"Generated prosite_scan_summary.tsv! PROSITE motifs found:\n{prosite_motifs}")

    #print("Running BLAST on protein sequences...")
    #blast_db = make_blast_db(out_dir, out_fasta)
    #blast_ref = select_blast_ref(records, out_dir)
    #blast_output = run_blast_search(out_dir, blast_ref, blast_db)
    #print("Successfully ran BLAST")

    #out_temp_fasta = temp_fasta(records, out_dir)
    #run_clustalo(out_dir, out_temp_fasta)
    #run_plotcon(out_dir)

    #out_motifs = run_prosite_scan(out_dir, out_fetch)
    #parse_prosite_output(out_dir, out_motifs)
    

if __name__ == "__main__":
    main()
