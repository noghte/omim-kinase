import json

n = 4  # Number of characters before and after the alignment position for the matched property

def parse_fasta_file(fasta_file_path):
    uniprot_info = {}
    with open(fasta_file_path, 'r') as file:
        current_id = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('|')
                if len(parts) > 1:
                    current_id = parts[0][1:] #parts[1]
                    sequence = next(file).strip()
                    sequence = sequence[sequence.index("{")+1:sequence.index("}")]
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
                        # elif char != '-':
                        else:
                            if in_flanking:
                                pos += 1
                            else:
                                if kinase_start == -1:
                                    kinase_start = pos
                                kinase_domain += char
                                pos += 1
                    if kinase_end == -1:  # If no flanking regions or at the end
                        kinase_end = pos - 1

                    # Create kinase_domain_alignment sequence
                    kinase_domain_alignment = ''.join(c for c in kinase_domain if c.isupper() or c == '-')

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
                        ],
                        "kinase_domain_alignment": {
                            "sequence": kinase_domain_alignment
                        }
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
                uniprot_id, from_aa, full_sequence_pos, to_aa = parts
                full_sequence_pos = int(full_sequence_pos)
                if uniprot_id in uniprot_info:
                    sequence = uniprot_info[uniprot_id]["sequence"]
                    
                    # Step 1: Find marker position (counting only letters, ignoring '(', ')', '-')
                    letter_count = 0
                    marker_pos = -1
                    for i, c in enumerate(sequence):
                        if c not in '( ) -':
                            letter_count += 1
                        if letter_count == full_sequence_pos:
                            marker_pos = i  # Save the position where we stop
                            break
                    
                    # Step 2: Count only uppercase letters and dashes to determine alignment_pos
                    alignment_pos = 0
                    inside_parentheses = False
                    for i, c in enumerate(sequence):
                        if c == '(':
                            inside_parentheses = True
                        elif c == ')':
                            inside_parentheses = False
                        elif not inside_parentheses and (c.isupper() or c == '-'):
                            alignment_pos += 1
                        
                        if i == marker_pos:
                            break
                    
                    # Determine location (kinase_domain or flanking_region)
                    location = "kinase_domain"
                    for flanking_region in uniprot_info[uniprot_id]["flanking_positions"]:
                        if flanking_region["start"] <= full_sequence_pos <= flanking_region["end"]:
                            location = "flanking_region"
                            alignment_pos = "Outside of the alignment"
                            break
                    
                    uniprot_info[uniprot_id]["substitutions"].append({
                        "full_sequence_pos": full_sequence_pos,
                        "alignment_pos": alignment_pos,
                        "from": from_aa,
                        "to": to_aa,
                        "location": location
                    })
    return uniprot_info

def write_json(output_file_path, uniprot_info):
    with open(output_file_path, 'w') as json_file:
        json.dump(list(uniprot_info.values()), json_file, indent=4)

def main():
    fasta_file_path = './kinsnps/subkinsnps.mma'
    # subs_file_path = './kinsnps/kinsnps_uid_subs_split.txt'
    subs_file_path_omim = './kinsnps/subkinsnps_uid_subs_split.txt'
    # subs_file_path_clinvar = './'
    output_file_path = 'kinsnps_allinfo.json'

    uniprot_info = parse_fasta_file(fasta_file_path)
    uniprot_info = parse_subs_file(subs_file_path_omim, uniprot_info)
    write_json(output_file_path, uniprot_info)

if __name__ == "__main__":
    main()
