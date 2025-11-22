#!.venv/bin/python3

import subprocess
import pandas as pd

def parse_clu_results(records, out_dir, out_clustalo, infoalign_results, out_prosite, size):
    infoalign_df = pd.read_csv(f"{out_dir}/{out_clustalo}/{infoalign_results}", sep = "\t", na_values = [""], header = None, names = ["Name", "% Change"])
    infoalign_df = infoalign_df.sort_values("% Change", ascending = False).reset_index(drop = True)
    top_df = infoalign_df.iloc[:size]
    top_df.to_csv(f"{out_dir}/{out_prosite}/top_{size}_from_clustalo.tsv", sep = "\t", header = True)
    list_top = top_df["Name"].tolist()
    top_records = list(filter(lambda x: x["uid"] in list_top, records))
    return top_records

def build_prosite_input(records_top, out_dir, out_prosite):
    prosite_input = []
    for entry in records_top:
        with open(f"{out_dir}/{out_prosite}/{entry['uid']}.fasta", "w") as file:
            file.write(f">{entry['uid']} {entry['description']} [{entry['species']}]\n")
            file.write(f"{entry['sequence']}\n\n")
        prosite_input.append(f"{entry['uid']}")
    return prosite_input

def run_prosite_scan(out_dir, out_prosite, prosite_input_list):
    for uid in prosite_input_list:
        prosite_input = f"{uid}.fasta"
        prosite_output = f"{uid}.patmatmotifs"
        subprocess.call(f"patmatmotifs -sequence {out_dir}/{out_prosite}/{prosite_input} -outfile {out_dir}/{out_prosite}/{prosite_output} -rformat excel -auto", shell = True, stdout = subprocess.DEVNULL)

def parse_prosite_output(out_dir, out_prosite, prosite_output_list):
    df_list = []
    for file in prosite_output_list:
        df = pd.read_csv(f"{out_dir}/{out_prosite}/{file}", sep = "\t", comment = "#", na_values = [""])
        df.columns = df.columns.str.strip()
        df_list.append(df)
    summary_df = pd.concat(df_list, join = "outer", ignore_index = True)
    summary_df.to_csv(f"{out_dir}/{out_prosite}/prosite_scan_summary.tsv", sep = "\t", header = True, index = False)
    prosite_motifs = list(set(summary_df["Motif"]))
    return summary_df, prosite_motifs

