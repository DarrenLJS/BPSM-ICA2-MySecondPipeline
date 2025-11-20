#!.venv/bin/python3

import subprocess

def build_clustalo_input(records, out_dir):
    records_sorted = list(sorted(records, key = lambda x: x["length"], reverse = True))
    if len(records_sorted) > 1000:
        records_sorted = records_sorted[:1000]
    clustalo_input = f"{out_dir}_clustalo_in.fasta"
    with open(f"{out_dir}/{clustalo_input}", "w") as file:
        for entry in records_sorted:
            file.write(f">{entry['uid']} {entry['description']} [{entry['species']}]\n")
            file.write(f"{entry['sequence']}\n\n")
    return clustalo_input

def run_clustalo(out_dir, clustalo_input_filename):
    clustalo_output = f"{out_dir}_clustalo_out.clu"
    subprocess.call(f"clustalo -i {out_dir}/{clustalo_input_filename} -o {out_dir}/{clustalo_output} --threads=4 --outfmt=clu", shell = True)
    return clustalo_output

def run_plotcon(out_dir, clustalo_output_filename):
    plotcon_graph = f"{out_dir}_plotcon_graph"
    subprocess.call(f"plotcon -sequence {out_dir}/{clustalo_output_filename} -graph png -goutfile {out_dir}/{plotcon_graph} -auto", shell = True)
    return f"{plotcon_graph}.1.png"

