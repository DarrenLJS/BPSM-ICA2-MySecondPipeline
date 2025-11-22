#!.venv/bin/python3

import requests
import sys

# run_esearch() takes in 3 args: protein family name, taxon ID and email.
# Function: Query NCBI Protein database using Entrez esearch REST API to get a list of unique IDs.
# If query failure or no results, print error message and return an empty list.
# Else, return a list of unique IDs.
def run_esearch(protein_family, taxon_id, email):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    query = f"({protein_family}) AND txid{taxon_id}[organism] NOT partial"
    params_esearch = {
        "db" : "protein", 
        "term" : query, 
        "retmax" : 5000, 
        "retmode" : "json", 
        "email" : email
    }
    response_esearch = requests.get(f"{base_url}/esearch.fcgi", params = params_esearch)
    if not response_esearch.ok:
        print(f"Error querying NCBI Protein! Error: {response_esearch.status_code}")
        return []

    search_result = response_esearch.json()
    ids_list = search_result.get("esearchresult", {}).get("idlist", [])
    if not ids_list:
        print(f"No result for {query} on NCBI Protein!")
        return []
    
    return ids_list

# run_efetch() takes in 3 args: list of unique IDs, out_dir and email.
# Function: With the list of unique IDs obtained from esearch, query NCBI Protein database using Entrez efetch REST API to get a raw FASTA file of all protein sequences.
# If query failure, print error message and exit system.
# Else, return the downloaded raw FASTA filename.
def run_efetch(ids_list, out_dir, email):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    ids_str = ",".join(ids_list)
    params_efetch = {
        "db" : "protein", 
        "id" : ids_str, 
        "rettype" : "fasta", 
        "retmode" : "text", 
        "email" : email
    }
    response_efetch = requests.post(f"{base_url}/efetch.fcgi", data = params_efetch)
    if not response_efetch.ok:
        print(f"Error fetching from NCBI Protein! Error: {response_efetch.status_code}")
        sys.exit("Please try again.")

    else:
        out_filename = f"{out_dir}_rawsequences.fasta"
        with open(f"./{out_dir}/{out_filename}", "w") as file:
            file.write(response_efetch.text)
        return out_filename

