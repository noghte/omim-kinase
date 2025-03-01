import json
from typing import Dict, List
import csv
import os
from typing import Dict, List, Optional
import re
protein_change_pattern = re.compile(r'^[A-Z]\d+[A-Z]$')

n = 4  # Number of characters before and after the alignment position for the matched property

def extract_sequence_info(sequence: str) -> tuple:
    """Extract kinase domain and flanking positions from sequence."""
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
        else:
            if in_flanking:
                pos += 1
            else:
                if kinase_start == -1:
                    kinase_start = pos
                kinase_domain += char
                pos += 1
                
    if kinase_end == -1:
        kinase_end = pos - 1
        
    return flanking_positions, kinase_domain, kinase_start, kinase_end

def parse_fasta_file(fasta_file_path: str) -> Dict:
    """Parse FASTA file and extract protein information."""
    uniprot_info = {}
    try:
        with open(fasta_file_path, 'r') as file:
            current_id = None
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    parts = line.split('|')
                    if len(parts) > 1:
                        current_id = parts[0][1:]
                        sequence = next(file).strip()
                        sequence = sequence[sequence.index("{")+1:sequence.index("}")]
                        
                        flanking_positions, kinase_domain, kinase_start, kinase_end = extract_sequence_info(sequence)
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
    except FileNotFoundError:
        print(f"Error: Could not find file {fasta_file_path}")
        return {}
    return uniprot_info

def get_alignment_position(sequence: str, marker_pos: int) -> int:
    """Calculate alignment position in sequence."""
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
    return alignment_pos

def parse_subs_file(subs_file_path: str, uniprot_info: Dict) -> Dict:
    """Parse substitutions file and update protein information."""
    try:
        with open(subs_file_path, 'r') as file:
            next(file)  # Skip header
            for line in file.read().splitlines():
                if not line:
                    continue
                    
                parts = line.split()
                if len(parts) == 4:
                    uniprot_id, from_aa, full_sequence_pos, to_aa = parts
                    full_sequence_pos = int(full_sequence_pos)
                    
                    if uniprot_id not in uniprot_info:
                        continue

                    sequence = uniprot_info[uniprot_id]["sequence"]
                    # Calculate letter_count first
                    current_count = 0
                    for i, c in enumerate(sequence):
                        if c not in '() -':
                            current_count += 1
                            if current_count == full_sequence_pos:
                                marker_pos = i
                                break
                    else:
                        marker_pos = -1
                    
                    alignment_pos = get_alignment_position(sequence, marker_pos)
                    
                    location = "kinase_domain"
                    for region in uniprot_info[uniprot_id]["flanking_positions"]:
                        if region["start"] <= full_sequence_pos <= region["end"]:
                            location = "flanking_region"
                            alignment_pos = "Outside of the alignment"
                            break
                    
                    uniprot_info[uniprot_id]["substitutions"].append({
                        "full_sequence_pos": full_sequence_pos,
                        "alignment_pos": alignment_pos,
                        "from": from_aa,
                        "to": to_aa,
                        "location": location,
                        "database": "OMIM"
                    })
    except FileNotFoundError:
        print(f"Error: Could not find file {subs_file_path}")
    return uniprot_info

def create_omim_uniprot_mapping() -> Dict[str, str]:
    """Create mapping between OMIM IDs and Uniprot IDs."""
    mapping = {}
    
    # Process found mappings
    with open('./data/omim_ids_found_with_uniprot.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            mapping[row['mimNumber']] = row['uniprot_id']
    
    # Process not found mappings
    with open('./data/omim_ids_notfound_with_uniprot.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if 'uniprot_id' in row and row['uniprot_id']:  # Only if uniprot_id exists
                mapping[row['mimNumber']] = row['uniprot_id']
    
    return mapping

def parse_clinvar_file(clinvar_file_path: str, uniprot_id: str, uniprot_info: Dict, omim_id: str) -> None:
    """Parse ClinVar file and add substitutions to uniprot_info."""
    print(f"Processing: Uniprot ID: {uniprot_id}, OMIM ID: {omim_id}")
    
    unmatched_file = './data/unmatched_protein_changes.csv'
    # Create/open the CSV file for unmatched patterns
    if not os.path.exists(unmatched_file):
        with open(unmatched_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['uniprot_id', 'omim_id', 'protein_change'])

    try:
        with open(clinvar_file_path, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                protein_change = row['Protein change']
                if not protein_change:  # Skip empty protein changes
                    continue
                if not protein_change_pattern.match(protein_change):
                    # Log only non-empty unmatched protein changes
                    with open(unmatched_file, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([uniprot_id, omim_id, protein_change])
                    continue
                try:
                    full_sequence_pos = int(''.join(c for c in protein_change if c.isdigit()))
                    from_aa = protein_change[0]  # First letter
                    to_aa = protein_change[-1]   # Last letter
                    
                    sequence = uniprot_info[uniprot_id]["sequence"]
                    # Calculate marker position
                    current_count = 0
                    for i, c in enumerate(sequence):
                        if c not in '() -':
                            current_count += 1
                            if current_count == full_sequence_pos:
                                marker_pos = i
                                break
                    else:
                        marker_pos = -1
                    
                    alignment_pos = get_alignment_position(sequence, marker_pos)
                    
                    # Determine location
                    location = "kinase_domain"
                    for region in uniprot_info[uniprot_id]["flanking_positions"]:
                        if region["start"] <= full_sequence_pos <= region["end"]:
                            location = "flanking_region"
                            alignment_pos = "Outside of the alignment"
                            break
                    
                    uniprot_info[uniprot_id]["substitutions"].append({
                        "full_sequence_pos": full_sequence_pos,
                        "alignment_pos": alignment_pos,
                        "from": from_aa,
                        "to": to_aa,
                        "location": location,
                        "database": "ClinVar"
                    })
                except (ValueError, KeyError) as e:
                    print(f"Error processing row in {clinvar_file_path}: {e}")
                    continue
    except FileNotFoundError:
        print(f"Warning: Could not find ClinVar file {clinvar_file_path}")

def add_clinvar_substitutions(uniprot_info: Dict) -> Dict:
    """Add ClinVar substitutions to uniprot_info."""
    omim_uniprot_mapping = create_omim_uniprot_mapping()
    
    # Create reverse mapping (Uniprot -> OMIM)
    uniprot_omim_mapping = {v: k for k, v in omim_uniprot_mapping.items()}
    
    for uniprot_id in uniprot_info:
        if uniprot_id in uniprot_omim_mapping:
            omim_id = uniprot_omim_mapping[uniprot_id]
            clinvar_file_path = f'./data/clinvar/{omim_id}.txt'
            parse_clinvar_file(clinvar_file_path, uniprot_id, uniprot_info, omim_id)
    
    return uniprot_info

if __name__ == "__main__":
    fasta_file_path = './kinsnps/subkinsnps.mma'
    subs_file_path_omim = './kinsnps/subkinsnps_uid_subs_split.txt'
    output_file_path = 'kinsnps_allinfo_twodbs.json'

    uniprot_info = parse_fasta_file(fasta_file_path)
    uniprot_info = parse_subs_file(subs_file_path_omim, uniprot_info)
    uniprot_info = add_clinvar_substitutions(uniprot_info)  # Add this line
    
    with open(output_file_path, 'w') as json_file:
        json.dump(list(uniprot_info.values()), json_file, indent=4)

    print("Done!")
