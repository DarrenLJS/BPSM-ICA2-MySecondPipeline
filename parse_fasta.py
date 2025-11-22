#!.venv/bin/python3

# parse_fasta() takes in 2 args: out_dir and raw FASTA filename.
# Function: Parse the downloaded raw FASTA file with all sequences into JSON format.
# Output: All protein sequences in JSON format.
def parse_fasta(out_dir, fasta_filename):
    with open(f"{out_dir}/{fasta_filename}", "r") as file:
        fasta_text = file.read()
    fasta_list = [entry.strip() for entry in fasta_text.split(">") if entry.strip()]
    
    def mapper(entry):
        entry = entry.replace("\n", "")
        parts = entry.split("]", 1)
        if len(parts) != 2:
            return None

        header, sequence = parts
        sequence = sequence.strip()
        header_parts = header.split("[")
        if len(header_parts) != 2:
            return None
        
        header_info, species = header_parts
        species = species.strip()
        info_parts = header_info.strip().split()
        if not info_parts:
            return None
        
        uid = info_parts[0]
        description = " ".join(info_parts[1:])
        return {
            "uid": uid, 
            "species": species, 
            "description": description, 
            "sequence": sequence, 
            "length": len(sequence)
        }
    
    records = list(map(lambda x: mapper(x), fasta_list))
    records = [record for record in records if record]
    return records
