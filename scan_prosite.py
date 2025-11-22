#!.venv/bin/python3

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# parse_clu_results() takes in 6 args: Parsed sequences to JSON, out_dir, out_clustalo, infoalign results file, out_prosite and user's input for scan size.
# Function: As per user's input, filter and select the top few sequences with high % similarity to the consensus sequence from infoalign results.
# Output: The top sequences in JSON and infoalign results.
def parse_clu_results(records, out_dir, out_clustalo, infoalign_results, out_prosite, size):
    infoalign_df = pd.read_csv(f"{out_dir}/{out_clustalo}/{infoalign_results}", sep = "\t", na_values = [""], header = None, names = ["Name", "% Change"])
    infoalign_df = infoalign_df.sort_values("% Change", ascending = False).reset_index(drop = True)
    top_df = infoalign_df.iloc[:size]
    top_df.to_csv(f"{out_dir}/{out_prosite}/top_{size}_from_clustalo.tsv", sep = "\t", header = True, index = False)
    list_top = top_df["Name"].tolist()
    top_records = list(filter(lambda x: x["uid"] in list_top, records))
    return top_records, top_df

# build_prosite_input() takes in 3 args: Top sequences in JSON, out_dir, out_prosite.
# Function: Prepare input FASTA files for patmatmotifs.
# Output: List of input FASTA files(filepaths) for patmatmotifs.
def build_prosite_input(records_top, out_dir, out_prosite):
    prosite_input = []
    for entry in records_top:
        with open(f"{out_dir}/{out_prosite}/{entry['uid']}.fasta", "w") as file:
            file.write(f">{entry['uid']} {entry['description']} [{entry['species']}]\n")
            file.write(f"{entry['sequence']}\n\n")
        prosite_input.append(f"{entry['uid']}")
    return prosite_input

# run_prosite_scan() takes in 3 args: out_dir, out_prosite and list of input FASTA files(filepaths) for patmatmotifs.
# Function: Run patmatmotifs command line for each input FASTA file, to scan for PROSITE motifs.
# Output: List of output patmatmotifs files.
def run_prosite_scan(out_dir, out_prosite, prosite_input_list):
    for uid in prosite_input_list:
        prosite_input = f"{uid}.fasta"
        prosite_output = f"{uid}.patmatmotifs"
        subprocess.call(f"patmatmotifs -sequence {out_dir}/{out_prosite}/{prosite_input} -outfile {out_dir}/{out_prosite}/{prosite_output} -rformat excel -auto", shell = True, stdout = subprocess.DEVNULL)

# parse_prosite_output() takes in 3 args: out_dir, out_prosite and list of output patmatmotifs files.
# Function: Parse output patmatmotifs files into a summary table, and plot a count plot of PROSITE motifs found.
# Output: Summary table of PROSITE motif locations per sequence, count plot of PROSITE motifs found, and list of unique PROSITE motifs found.
def parse_prosite_output(out_dir, out_prosite, prosite_output_list):
    df_list = []
    for file in prosite_output_list:
        df = pd.read_csv(f"{out_dir}/{out_prosite}/{file}", sep = "\t", comment = "#", na_values = [""])
        df.columns = df.columns.str.strip()
        df_list.append(df)
    summary_df = pd.concat(df_list, join = "outer", ignore_index = True)
    summary_df.to_csv(f"{out_dir}/{out_prosite}/prosite_scan_summary.tsv", sep = "\t", header = True, index = False)
    
    motif_counts = summary_df["Motif"].value_counts()
    plt.figure(figsize = (12, 6))
    ax = sns.countplot(x = "Motif", data = summary_df, hue = "Motif", order = motif_counts.index)
    plt.title("Count Plot of PROSITE Motifs Found")
    for container in ax.containers:
        ax.bar_label(container, fmt = "%d")
    plt.tight_layout()
    
    print("Generating count_plot_motifs.png...")
    plt.savefig(f"{out_dir}/{out_prosite}/count_plot_motifs.png")
    print("Outputting Count Plot of PROSITE Motifs Found to terminal. Check plot, and close to continue!")
    plt.show()
    plt.close()
    
    temp_df = summary_df.copy()
    temp_df["Start, End, Score"] = temp_df.apply(lambda x: [x["Start"], x["End"], ["Score"]], axis = 1)
    prosite_df = temp_df.pivot_table(index = "SeqName", columns = "Motif", values = "Start, End, Score", aggfunc = list, fill_value = "x")
    prosite_df.to_csv(f"{out_dir}/{out_prosite}/prosite_locations.tsv", sep = "\t", header = True, index = False)
    return summary_df, motif_counts.index.tolist()

