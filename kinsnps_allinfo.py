import json

def parse_fasta_file(fasta_file_path):
    uniprot_info = {}
    with open(fasta_file_path, 'r') as file:
        current_id = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('|')
                if len(parts) > 1:
                    current_id = parts[1]
                    sequence = next(file).strip()
                    sequence = sequence[sequence.index("{")+1:sequence.index("}")]
                    clean_sequence = sequence.replace('(', '').replace(')', '').replace('-', '')
                    flanking_positions = []
                    kinase_domain = ""
                    kinase_start = -1
                    kinase_end = -1
                    pos = 1
                    in_flanking = False
                    start_pos = -1
                    for i, char in enumerate(sequence):
                        if char == '(':
                            in_flanking = True
                            start_pos = pos
                            if kinase_domain:
                                kinase_end = pos - 1
                        elif char == ')':
                            in_flanking = False
                            flanking_positions.append({"start": start_pos, "end": pos-1})
                            if kinase_start == -1:
                                kinase_start = pos
                        elif char != '-':
                            if in_flanking:
                                pos += 1
                            else:
                                if kinase_start == -1:
                                    kinase_start = pos
                                kinase_domain += char
                                pos += 1
                    if kinase_end == -1:  # If no flanking regions or at the end
                        kinase_end = pos - 1
                    uniprot_info[current_id] = {
                        "uniprot_id": current_id,
                        "sequence": sequence,
                        "substitutions": [],
                        "flanking_positions": flanking_positions,
                        "kinase_domain": {
                            "sequence": kinase_domain.replace("-", ""),
                            "start": kinase_start,
                            "end": kinase_end
                        },
                        "kinase_motifs": [
                            {"name": "", "start": -1, "end": -1}
                        ]
                    }
    return uniprot_info

def parse_subs_file(subs_file_path, uniprot_info):
    with open(subs_file_path, 'r') as file:
        next(file)  # Skip header
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) == 4:
                uniprot_id, from_aa, pos, to_aa = parts
                pos = int(pos)
                if uniprot_id in uniprot_info:
                    uniprot_info[uniprot_id]["substitutions"].append({
                        "pos": pos,
                        "from": from_aa,
                        "to": to_aa
                    })
    return uniprot_info

def write_json(output_file_path, uniprot_info):
    with open(output_file_path, 'w') as json_file:
        json.dump(list(uniprot_info.values()), json_file, indent=4)

def main():
    fasta_file_path = './kinsnps/kinsnps_subs.mma'
    subs_file_path = './kinsnps/kinsnps_uid_subs_split.txt'
    output_file_path = 'kinsnps_allinfo.json'

    uniprot_info = parse_fasta_file(fasta_file_path)
    uniprot_info = parse_subs_file(subs_file_path, uniprot_info)
    write_json(output_file_path, uniprot_info)

if __name__ == "__main__":
    main()
