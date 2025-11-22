#!.venv/bin/python3

import subprocess

# build_clustalo_input() takes in 4 args: Parsed sequences to JSON, out_dir, out_clustalo and user's input for sample size.
# Function: Prepare an input FASTA file for ClustalO.
# As per user's input size, sort the sequences by length in descending order, and select top few longest sequences for conservation analysis.
# Output: Input FASTA filename for Clustalo.
def build_clustalo_input(records, out_dir, out_clustalo, size):
    records_sorted = list(sorted(records, key = lambda x: x["length"], reverse = True))
    if len(records_sorted) > size:
        records_sorted = records_sorted[:size]
    clustalo_input = f"{out_clustalo}_in.fasta"
    with open(f"{out_dir}/{out_clustalo}/{clustalo_input}", "w") as file:
        for entry in records_sorted:
            file.write(f">{entry['uid']} {entry['description']} [{entry['species']}]\n")
            file.write(f"{entry['sequence']}\n\n")
    return clustalo_input

# run_clustalo() takes in 3 args: out_dir, out_clustalo and input FASTA filename for ClustalO.
# Function: Run clustalo command line.
# Output: Multi-Sequence Alignment(MSA) output file with .clu.
def run_clustalo(out_dir, out_clustalo, clustalo_input_filename):
    clustalo_output = f"{out_clustalo}_out.clu"
    subprocess.call(f"clustalo -i {out_dir}/{out_clustalo}/{clustalo_input_filename} -o {out_dir}/{out_clustalo}/{clustalo_output} --threads=4 --outfmt=clu", shell = True, stdout = subprocess.DEVNULL)
    return clustalo_output

# run_plotcon() takes in 3 args: out_dir, out_clustalo and MSA output.
# Function: Run plotcon command lines to generate Plotcon graph to terminal and save graph as .png file.
# Run plotcon command line to generate Plotcon raw data and save as .txt file.
# Output: Plotcon graph and raw data files.
def run_plotcon(out_dir, out_clustalo, clustalo_output_filename):
    plotcon_graph = f"{out_clustalo}_plotcon_graph"
    plotcon_data = f"{out_clustalo}_plotcon_data"
    print(f"Outputting {plotcon_graph} to terminal. Check graph, and close to continue!")
    subprocess.call(f"plotcon -sequence {out_dir}/{out_clustalo}/{clustalo_output_filename} -graph x11 -auto", shell = True, stdout = subprocess.DEVNULL)
    subprocess.call(f"plotcon -sequence {out_dir}/{out_clustalo}/{clustalo_output_filename} -graph png -goutfile {out_dir}/{out_clustalo}/{plotcon_graph} -auto", shell = True, stdout = subprocess.DEVNULL)
    subprocess.call(f"plotcon -sequence {out_dir}/{out_clustalo}/{clustalo_output_filename} -graph data -goutfile {out_dir}/{out_clustalo}/{plotcon_data} -auto", shell = True, stdout = subprocess.DEVNULL)
    return f"{plotcon_graph}.1.png", f"{plotcon_data}1.dat"

# get_consensus() takes in 3 args: out_dir, out_clustalo and MSA output.
# Function: Run cons command line to obtain consensus sequence from MSA output
# Output: Consensus sequence stored in FASTA format.
def get_consensus(out_dir, out_clustalo, clustalo_output_filename):
    consensus = f"{out_clustalo}_consensus.fasta"
    subprocess.call(f"cons -sequence {out_dir}/{out_clustalo}/{clustalo_output_filename} -outseq {out_dir}/{out_clustalo}/{consensus} -auto", shell = True, stdout = subprocess.DEVNULL)
    return consensus

# get_infoalign() takes in 3 args: out_dir, out_clustalo and MSA output.
# Function: Run infoalign command line to obtain % similarity of each sequence to consensus sequence.
# Output: infoalign results .txt file.
def get_infoalign(out_dir, out_clustalo, clustalo_output_filename):
    infoalign_results = f"{out_clustalo}_infoalign.txt"
    subprocess.call(f"infoalign {out_dir}/{out_clustalo}/{clustalo_output_filename} -outfile {out_dir}/{out_clustalo}/{infoalign_results} -only -name -change -auto", shell = True, stdout = subprocess.DEVNULL)
    return infoalign_results

