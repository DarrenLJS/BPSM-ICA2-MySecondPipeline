#!.venv/bin/python3

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# make_blast_db() takes in 3 args: out_dir, raw FASTA file of all protein sequences and blast_db.
# Function: Run makeblastdb command line to construct a BLAST database.
# Output: Custom BLAST database name.
def make_blast_db(out_dir, fasta_filename, blast_db):
    subprocess.call(f"makeblastdb -in {out_dir}/{fasta_filename} -dbtype prot -out {out_dir}/{blast_db}/{blast_db}", shell = True)

# select_blast_ref() takes in 4 args: Top sequences in JSON, top sequences in infoalign results, out_dir and blast_db.
# Function: Select the sequence with the highest % similarity with the consensus sequence and save the sequence in a FASTA file.
# Output: BLAST reference sequence FASTA file.
def select_blast_ref(top_records, top_df, out_dir, blast_db):
    ref_seq = top_df.iloc[0]["Name"]
    blast_ref = list(filter(lambda x: x["uid"] == ref_seq, top_records))[0]
    ref_file = f"{blast_ref['uid']}_blast_ref.fasta"
    with open(f"{out_dir}/{blast_db}/{ref_file}", "w") as file:
        file.write(f">{blast_ref['uid']} {blast_ref['description']} [{blast_ref['species']}]\n")
        file.write(f"{blast_ref['sequence']}\n\n")
    return ref_file

# run_blast_search takes in 3 args: out_dir, BLAST reference sequence FASTA file and blast_db.
# Function: Run blastp command line to conduct BLAST analysis of protein database with a protein sequence query.
# Output: BLAST output file.
def run_blast_search(out_dir, blastref_filename, blast_db):
    blast_output = f"{out_dir}_blastoutput.out"
    subprocess.call(f"blastp -db {out_dir}/{blast_db}/{blast_db} -query {out_dir}/{blast_db}/{blastref_filename} -outfmt 7 > {out_dir}/{blast_db}/{blast_output}", shell = True)
    return blast_output

# arse_blast_output() takes in 3 args: out_dir, blast_db and BLAST output file.
# Function: Parse BLAST output file to a pandas dataframe, and visualise BLAST output by plotting a scatterplot of BLAST hits (Alignment Length against % Identity).
# Output: Scatterplot of BLAST hits and parsed BLAST output to pandas dataframe.
def parse_blast_output(out_dir, blast_db, blastresults_filename):
    fields = ["query", "subject", "pid", "length", "mismatch", "gap_open", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
    blast_df = pd.read_csv(f"{out_dir}/{blast_db}/{blastresults_filename}", sep = "\t", comment = "#", names = fields)

    plt.figure(figsize = (12, 6))
    sns.scatterplot(x = "pid", y = "length", data = blast_df, hue = "bitscore", palette = "Blues", s = 40, edgecolor = "black")
    plt.title("Scatterplot of BLAST Hits")
    plt.xlabel("% Identity")
    plt.ylabel("Alignment Length")
    plt.grid(color = "lightgray", linestyle = "--")
    plt.tight_layout()
    
    print("Generating scatterplot_blast.png...")
    plt.savefig(f"{out_dir}/{blast_db}/scatterplot_blast.png")
    print("Outputting Scatterplot of BLAST Hits to terminal. Check plot, and close to continue!")
    plt.show()
    plt.close()
    return blast_df

