#!.venv/bin/python3

import subprocess
import requests

def validate_protein(protein_family):
    base_url = "https://www.ebi.ac.uk/ebisearch/ws/rest/pfam"
    params = {
        "query" : protein_family
    }
    headers = {"Accept" : "application/json"}
    response = requests.get(base_url, params = params, headers = headers)
    if not response.ok:
        raise ValueError("Error querying EBI Search! Error: {response.status_code}")

    data = response.json()
    results = data.get("entries", [])
    if not results:
        raise ValueError(f"No result for {protein_family} on EBI Search!")

    pfam_id = results[0].get("acc")
    pfam_name = results[0].get("id")
    return pfam_id, pfam_name

def validate_taxon(taxon_group, email):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    params_esearch = {
        "db" : "taxonomy", 
        "term" : taxon_group, 
        "retmax" : 20, 
        "retmode" : "json", 
        "email" : email
    }
    response_esearch = requests.get(f"{base_url}/esearch.fcgi", params = params_esearch)
    if not response_esearch.ok:
        raise ValueError("Error querying NCBI Taxonomy through esearch! Error: {response_esearch.status_code}")

    search_result = response_esearch.json()
    id_list = search_result.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        raise ValueError(f"No result for {taxon_group} on NCBI Taxonomy!")
    
    taxon_id = id_list[0]
    params_esummary = {
        "db" : "taxonomy", 
        "id" : taxon_id, 
        "retmode" : "json", 
        "email" : email
    }
    response_esummary = requests.get(f"{base_url}/esummary.fcgi", params = params_esummary)
    if not response_esummary.ok:
        raise ValueError("Error querying NCBI Taxonomy through esummary! Error: {response_esummary.status_code}")

    summary = response_esummary.json()
    taxon_name = summary.get("result", {}).get(taxon_id, {}).get("scientificname")
    if not taxon_name:
        raise ValueError(f"No result for {taxon_id} on NCBI Taxonomy!")

    return taxon_id, taxon_name

