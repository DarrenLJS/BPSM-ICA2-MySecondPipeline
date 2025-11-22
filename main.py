#!.venv/bin/python3

import glob
import json
import subprocess
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from input_handler import validate_protein, validate_taxon
from fetch_sequence import run_esearch, run_efetch
from parse_fasta import parse_fasta
from conservation_analysis import build_clustalo_input, run_clustalo, run_plotcon, get_consensus, get_infoalign
from scan_prosite import parse_clu_results, build_prosite_input, run_prosite_scan, parse_prosite_output
from blast_analysis import make_blast_db, select_blast_ref, run_blast_search, parse_blast_output

# main() runs the pipeline.
def main():
    
    # This while loop prompts the user for an email. If the email is invalid, reprompt the user, until a valid email is inputted.
    while True:
        email = input("Your email (needed for NCBI URL request): ").strip()
        if re.search(
            r"^[a-zA-Z0-9.!#$%&'*+\/=?^_`{|}~-]+@[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?(?:\.[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?)*$", 
            email, re.IGNORECASE
        ):
            break
        else:
            print("\nPlease input a valid email.")
    
    # This while loop prompts the user for a protein family, which is to be validated by validate_protein().
    # If input is invalid, raise an error and reprompt user, until a valid protein family is inputted.
    while True:
        try:
            pfam_input = input("Protein family: ").strip()
            pfam_id, pfam_name = validate_protein(pfam_input)
            break
        except Exception as e:
            print(f"\nPlease input a valid protein family. {e}")
    
    # This while loop prompts user for a taxonomic group, which is to be validated by validate_taxon().
    # If input is invalid, raise an error and reprompt user, until a valid taxon group is inputted.
    while True:
        try:
            taxon_input = input("Taxon group: ").strip()
            taxon_id, taxon_name = validate_taxon(taxon_input, email)
            break
        except Exception as e:
            print(f"\nPlease input a valid taxon group. {e}")
    
    # When all inputs are valid, print these statements.
    print(f"\nSuccess! Your protein family: {pfam_name}\tPfam ID: {pfam_id}")
    print(f"Success! Your taxon group: {taxon_name}\tTaxon ID: {taxon_id}\n")
    
    # Query NCBI Protein database using Entrez esearch to get a list of unique IDs, through run_esearch().
    # run_esearch() takes in 3 args: protein family name, taxon ID and email. Output: List of unique IDs.
    # If list of unique IDs is empty, exit system with "no result" message.
    print("Searching NCBI Protein...")
    ids_list = run_esearch(pfam_name, taxon_id, email)
    if not ids_list:
        sys.exit("\nNo UIDs could be retrieved from NCBI Protein! Please try again.")
    
    # If there are more than 1000 entries in the dataset, prompt user for confirmation.
    # This while loop prompts user for yes/no. If invalid input, raise an error and reprompt user, until valid.
    # If no, exit system. Else, continue.
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
    
    # Make an output directory named after protein family and taxon group. Store directory name as out_dir.
    print("Protein UIDs successfully retrieved!\n")
    out_dir = f"{pfam_name}_{taxon_name}"
    subprocess.call(f"rm -rf {out_dir}", shell = True)
    subprocess.call(f"mkdir -p {out_dir}", shell = True)
    
    # With the list of unique IDs obtained from esearch, query NCBI Protein database using Entrez efetch to get a raw FASTA file of all protein sequences, through run_efetch().
    # run_efetch() takes in 3 args: list of unique IDs, out_dir and email. Output: Raw FASTA filename.
    print("Fetching all protein sequences in FASTA format...")
    out_fasta = run_efetch(ids_list, out_dir, email)
    print(f"FASTA sequences successfully downloaded to {out_fasta}\n")
    
    # Parse raw FASTA file into JSON format, through parse_fasta().
    # parse_fasta() takes in 2 args: out_dir and raw FASTA filename. Output: All protein sequences in JSON format.
    print(f"Parsing {out_fasta}...")
    records = parse_fasta(out_dir, out_fasta)
    with open(f"{out_dir}/{out_dir}_records.json", "w") as file:
        json.dump(records, file, indent = 2)
    print(f"Successfully parsed {out_fasta}\n")
    
    # This while loop prompts user for an integer, for the number of sequences to be analysed with ClustalO.
    # Input needs to be an integer <= Total number of sequences. If invalid, raise an error and reprompt user, until valid.
    # If input > 2000 and exceeds recommended range, prompt user for a confirmation (yes/no). If invalid input, raise an error and reprompt user, until valid.
    # If no, reprompt user for another integer. Else, continue.
    while True:
        try:
            print(f"Please input a number for the number of sequences to be analysed with ClustalO. Maximum: {len(records)}")
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
    
    # Show user's input by printing these statements.
    print(f"\nYour sample size: {conservation_analysis_size}")
    print(f"The {conservation_analysis_size} longest protein sequences will be used for conservation analysis.\n")
    
    # Perform conservation analysis with ClustalO, through build_clustalo_input() and run_clustalo().
    # Make an output directory for ClustalO and store the directory name as out_clustalo.
    # build_clustalo_input() takes in 4 args: Parsed sequences to JSON, out_dir, out_clustalo and user's input for sample size. Output: Input FASTA filename for ClustalO.
    # run_clustalo() takes in 3 args: out_dir, out_clustalo and input FASTA filename for ClustalO. Output: Multi-Sequence Alignment(MSA) output filename with .clu.
    print("Performing conservation analysis with ClustalO and Plotcon...")
    out_clustalo = f"{out_dir}_clustalo"
    subprocess.call(f"mkdir -p {out_dir}/{out_clustalo}", shell = True)
    clustalo_input = build_clustalo_input(records, out_dir, out_clustalo, conservation_analysis_size)
    clustalo_output = run_clustalo(out_dir, out_clustalo, clustalo_input)
    print(f"ClustalO output successfully generated! Generated {clustalo_output}")

    # With the Multi-Sequence Alignment(MSA) output, visualise data with Plotcon, through run_plotcon().
    # run_plotcon() takes in 3 args: out_dir, out_clustalo and MSA output. Output: Plotcon graph and raw data filenames.
    # Obtain consensus sequence with cons, through get_consensus().
    # get_consensus() takes in 3 args: out_dir, out_clustalo and MSA output. Output: Consensus sequence FASTA filename.
    # Obtain % similarity of sequences to consensus sequence with infoalign, through get_infoalign().
    # get_infoalign() takes in 3 args: out_dir, out_clustalo and MSA output. Output: infoalign results .txt filename.
    plotcon_graph, plotcon_text = run_plotcon(out_dir, out_clustalo, clustalo_output)
    print(f"Plotcon output successfully generated! Generated {plotcon_graph}, {plotcon_text}")
    consensus = get_consensus(out_dir, out_clustalo, clustalo_output)
    print(f"Successfully obtained consensus sequence! Generated {consensus}")
    infoalign_results = get_infoalign(out_dir, out_clustalo, clustalo_output)
    print(f"Successfully obtained infoalign results for {clustalo_output}! Generated {infoalign_results}\n")
    
    # This while loop prompts user for an integer, for the number of sequences to be scanned for PROSITE motifs.
    # Input needs to be an integer <= Number of selected sequences for conservation analysis. If invalid, raise an error and reprompt user, until valid.
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
    
    # Show user's input by printing these statements.
    print(f"\nYour scan size: {scan_size}")
    print(f"The {scan_size} protein sequences with the highest similarity to the consensus sequence will be scanned for PROSITE motifs.\n")
    
    # Make an output directory for PROSITE scan and store the directory name as out_prosite.
    # Filter for the top sequences with high similarities to the consensus sequence, through parse_clu_results().
    # parse_clu_results() takes in 6 args: Parsed sequences to JSON, out_dir, out_clustalo, infoalign results file, out_prosite and user's input for scan size.
    # Output: The top sequences in JSON and infoalign results.
    out_prosite = f"{out_dir}_prosite"
    subprocess.call(f"mkdir -p {out_dir}/{out_prosite}", shell = True)
    top_records, top_df = parse_clu_results(records, out_dir, out_clustalo, infoalign_results, out_prosite, scan_size)
    
    # Scan protein sequences for PROSITE motifs using patmatmotifs, through build_prosite_input() and run_prosite_scan().
    # build_prosite_input() takes in 3 args: Top sequences in JSON, out_dir, out_prosite. Output: List of input FASTA files(filepaths) for patmatmotifs.
    # run_prosite_scan() takes in 3 args: out_dir, out_prosite and list of input FASTA files(filepaths) for patmatmotifs. Output: List of output patmatmotifs files.
    print("Scanning protein sequences for PROSITE motifs...")
    prosite_input = build_prosite_input(top_records, out_dir, out_prosite)
    run_prosite_scan(out_dir, out_prosite, prosite_input)
    prosite_output = list(map(lambda x: x + ".patmatmotifs", prosite_input))
    print(f"PROSITE output successfully generated! Generated {len(prosite_output)} .patmatmotifs files.")

    # Parse output patmatmotifs files into a summary table, and plot a count plot of PROSITE motifs found, through parse_prosite_output().
    # parse_prosite_output() takes in 3 args: out_dir, out_prosite and list of output patmatmotifs files.
    # Output: Summary table of PROSITE motif locations per sequence, count plot of PROSITE motifs found, and list of unique PROSITE motifs found.
    prosite_summary_df, prosite_motifs = parse_prosite_output(out_dir, out_prosite, prosite_output)
    print(f"Generated prosite_scan_summary.tsv! PROSITE motifs found:\n{prosite_motifs}\n")
    
    # Make an output directory for BLAST analysis and store the directory name as blast_db.
    # Construct a BLAST database using makeblastdb, through make_blast_db().
    # make_blast_db() takes in 3 args: out_dir, raw FASTA file of all protein sequences and blast_db. Output: Custom BLAST database name.
    print("Running BLAST on protein sequences...")
    blast_db = f"{out_dir}_blast"
    subprocess.call(f"mkdir -p {out_dir}/{blast_db}", shell = True)
    make_blast_db(out_dir, out_fasta, blast_db)

    # Select and construct a BLAST reference sequence FASTA file using the sequence with the highest similarity to the consensus sequence from ClustalO, through select_blast_ref().
    # select_blast_ref() takes in 4 args: Top sequences in JSON, top sequences in infoalign results, out_dir and blast_db. Output: BLAST reference sequence FASTA file.
    # Run BLAST analysis with blastp, through run_blast_search().
    # run_blast_search takes in 3 args: out_dir, BLAST reference sequence FASTA file and blast_db. Output: BLAST output file.
    blast_ref = select_blast_ref(top_records, top_df, out_dir, blast_db)
    print(f"Using {blast_ref} as BLAST reference sequence...")
    blast_output = run_blast_search(out_dir, blast_ref, blast_db)
    print(f"Successfully completed BLAST analysis! Generated {blast_output}")

    # Visualise BLAST output by plotting a scatterplot of BLAST hits (Alignment Length against % Identity), through parse_blast_output().
    # parse_blast_output() takes in 3 args: out_dir, blast_db and BLAST output file. Output: Scatterplot of BLAST hits and parsed BLAST output to pandas dataframe.
    blast_df = parse_blast_output(out_dir, blast_db, blast_output)

    print(f"\nFull analysis successfully completed! Check {out_dir} for all outputs.\n")
    

if __name__ == "__main__":
    main()
