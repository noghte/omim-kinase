# Validate the `subkinsnps_allinfo.json` file with the following steps:
# 1. From the Sequence: Delete paranthesis and dashes
# 2. Find the "from" character that is in the "full_seuqnce_pos" and see if it  corresponds with the "sequence"
# 3. Check if the "from" corresponds to  the "alignment_pos" in the "kinase_domain_alignment" 

import json

# Conversion table for amino acids from 3-letter to 1-letter codes
amino_acid_map = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "ASX": "B", "GLX": "Z"
}

def remove_parentheses_and_dashes(sequence):
    return sequence.replace("(", "").replace(")", "").replace("-", "")

def validate_substitutions(data):
    errors = []
    uniprot_errors = set()
    
    for entry in data:
        uniprot_id = entry['uniprot_id']
        # Step 1: Clean the sequence by removing parentheses and dashes
        cleaned_sequence = remove_parentheses_and_dashes(entry["sequence"])
        
        for substitution in entry["substitutions"]:
            full_sequence_pos = substitution["full_sequence_pos"] - 1  # Convert to 0-based index
            alignment_pos = substitution["alignment_pos"] - 1  # Convert to 0-based index
            from_residue = amino_acid_map.get(substitution["from"].upper())

            if not from_residue:
                errors.append({"UniprotId": uniprot_id, "Error": f"Invalid amino acid code {substitution['from']}"})
                uniprot_errors.add(uniprot_id)
                continue

            # Step 2: Validate the 'from' residue with the full sequence
            try:
                sequence_residue = cleaned_sequence[full_sequence_pos]
            except IndexError:
                errors.append({"UniprotId": uniprot_id, "Error": f"Index out of range in sequence at full_sequence_pos {substitution['full_sequence_pos']}"})
                uniprot_errors.add(uniprot_id)
                continue

            if sequence_residue.upper() != from_residue:
                errors.append({"UniprotId": uniprot_id, "Error": f"Mismatch in sequence at full_sequence_pos {substitution['full_sequence_pos']}: expected {from_residue}, found {sequence_residue}"})
                uniprot_errors.add(uniprot_id)

            # Step 3: Validate the 'from' residue with the kinase domain alignment if applicable
            if "kinase_domain" in entry:
                kinase_domain_residue = entry["kinase_domain_alignment"]["sequence"][alignment_pos]

                if kinase_domain_residue.upper() != from_residue:
                    errors.append({"UniprotId": uniprot_id, "Error": f"Mismatch in kinase domain alignment at alignment_pos {substitution['alignment_pos']}: expected {from_residue}, found {kinase_domain_residue}"})
                    uniprot_errors.add(uniprot_id)
    
    return errors, uniprot_errors

# Load the JSON file
with open('kinsnps_allinfo.json', 'r') as file:
    data = json.load(file)

# Perform validation
errors, uniprot_errors = validate_substitutions(data)

# Output results
if errors:
    for error in errors:
        print(error)
    print(f"\nTotal UniProts with errors: {len(uniprot_errors)} out of {len(data)}")
else:
    print("All substitutions are valid!")
    print(f"\nTotal UniProts checked: {len(data)}")