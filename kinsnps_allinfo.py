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
                    current_id = parts[1]
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
                            "sequence": kinase_domain,#.replace("-", ""),
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
                    # Determine location
                    location = "kinase_domain"
                    for flanking_region in uniprot_info[uniprot_id]["flanking_positions"]:
                        if flanking_region["start"] <= full_sequence_pos <= flanking_region["end"]:
                            location = "flanking_region"
                            alignment_pos = full_sequence_pos
                            break

                    if location == "kinase_domain":
                        # Calculate the length of all preceding flanking regions
                        preceding_flanking_length = sum(
                            (region["end"] - region["start"] + 1) 
                            for region in uniprot_info[uniprot_id]["flanking_positions"]
                            if region["end"] < full_sequence_pos
                        )

                        # Calculate alignment_pos
                        sequence = uniprot_info[uniprot_id]["sequence"]
                        upper_count = sum(1 for c in sequence[:full_sequence_pos] if c.isupper() or c == '(' or c == ')' or c == '-')
                        alignment_pos = upper_count - preceding_flanking_length

                    # Create the matched property
                    if location == "flanking_region":
                        sequence = uniprot_info[uniprot_id]["sequence"]
                        flanking_sequences = ''.join(c for c in sequence if c.isupper() or c == '(' or c == ')')
                        matched_start = max(0, alignment_pos - 1 - n)
                        matched_end = min(len(flanking_sequences), alignment_pos - 1 + n + 1)
                        matched_seq = flanking_sequences[matched_start:matched_end]
                    else:
                        cleaned_kinase_seq = uniprot_info[uniprot_id]["kinase_domain_alignment"]["sequence"]
                        matched_start = max(0, alignment_pos - 1 - n)
                        matched_end = min(len(cleaned_kinase_seq), alignment_pos - 1 + n + 1)
                        matched_seq = cleaned_kinase_seq[matched_start:matched_end]

                    uniprot_info[uniprot_id]["substitutions"].append({
                        "full_sequence_pos": full_sequence_pos,
                        "alignment_pos": alignment_pos,
                        "from": from_aa,
                        "to": to_aa,
                        "matched": matched_seq,
                        "location": location
                    })
    return uniprot_info

def write_json(output_file_path, uniprot_info):
    with open(output_file_path, 'w') as json_file:
        json.dump(list(uniprot_info.values()), json_file, indent=4)

def main():
    fasta_file_path = './kinsnps/kinsnps_subs.mma'
    # subs_file_path = './kinsnps/kinsnps_uid_subs_split.txt'
    subs_file_path = './kinsnps/subkinsnps_uid_subs_split.txt'
    output_file_path = 'kinsnps_allinfo.json'

    uniprot_info = parse_fasta_file(fasta_file_path)
    uniprot_info = parse_subs_file(subs_file_path, uniprot_info)
    write_json(output_file_path, uniprot_info)

if __name__ == "__main__":
    main()