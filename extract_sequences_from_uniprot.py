import requests
import csv
import time
import os

def get_protein_data(uniprot_id):
    time.sleep(0.1)
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    headers = {
        'Accept': 'application/json'
    }
    
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        return response.json()
    elif response.status_code == 404:
        print(f"UniProt ID {uniprot_id} not found.")
        return None
    else:
        print(f"Error {response.status_code} for UniProt ID {uniprot_id}.")
        return None

def get_isoform_sequences(uniprot_id):
    time.sleep(0.1)
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}/isoforms"
    headers = {
        'Accept': 'application/json'
    }
    
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        if isinstance(data, list) and len(data) > 0:
            data.sort(key=lambda x: int(x['accession'].split('-')[-1]) if '-' in x['accession'] else 0)
            return data
        else:
            return None
    else:
        return None

def sequence_matches_position(sequence, position, wt_amino_acid):
    try:
        return sequence[int(position) - 1] == wt_amino_acid
    except IndexError:
        return False

def read_tsv_and_generate_fasta(input_tsv, output_fasta, mode):
    processed_uniprot_ids = set()  # To track processed UniProt IDs
    
    with open(input_tsv, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        tsv_data = {}
        for row in reader:
            uniprot_id = row['Uniprot ID']
            if uniprot_id not in tsv_data:
                tsv_data[uniprot_id] = []
            tsv_data[uniprot_id].append(row)
        
        with open(output_fasta, mode) as fasta_file:
            for uniprot_id, positions in tsv_data.items():
                if uniprot_id in processed_uniprot_ids:
                    continue  # Skip if this UniProt ID has already been processed

                print(f"Processing Uniprot ID: {uniprot_id}")
                protein_data = get_protein_data(uniprot_id)
                
                if protein_data:
                    # Check the main protein sequence first
                    sequence = protein_data.get('sequence', {}).get('sequence', '')
                    for pos in positions:
                        if sequence_matches_position(sequence, pos['Position'], pos['WT Amino Acid']):
                            # Write the main protein sequence to the FASTA file
                            fasta_header = (
                                f">{protein_data['accession']}|{protein_data['id']}|{protein_data['protein']['recommendedName']['fullName']['value']}|"
                                f"GN={protein_data['gene'][0]['name']['value']}|OS={next((n['value'] for n in protein_data['organism']['names'] if n['type'] == 'scientific'), 'N/A')}|OX={protein_data['organism']['taxonomy']}"
                            )
                            fasta_file.write(f"{fasta_header}\n{sequence}\n\n")
                            processed_uniprot_ids.add(uniprot_id)
                            break  # Skip isoform checking if main sequence matches any position

                    # If the main sequence does not match any position, check isoforms
                    if uniprot_id not in processed_uniprot_ids:
                        isoforms = get_isoform_sequences(uniprot_id)
                        if isoforms:
                            for isoform in isoforms:
                                sequence = isoform.get('sequence', {}).get('sequence', '')
                                for pos in positions:
                                    if sequence_matches_position(sequence, pos['Position'], pos['WT Amino Acid']):
                                        # Write the isoform sequence to the FASTA file
                                        fasta_header = (
                                            f">{isoform['accession']}|{isoform['id']}|{isoform['protein']['recommendedName']['fullName']['value']}|"
                                            f"GN={isoform['gene'][0]['name']['value']}|OS={next((n['value'] for n in isoform['organism']['names'] if n['type'] == 'scientific'), 'N/A')}|OX={isoform['organism']['taxonomy']}"
                                        )
                                        fasta_file.write(f"{fasta_header}\n{sequence}\n\n")
                                        processed_uniprot_ids.add(uniprot_id)
                                        break  # Stop after the first matching isoform

def report_missings(tsv_file, fasta_file):
    tsv_uniprot_ids = set()

    # Step 1: Extract UniProt IDs from the TSV file
    with open(tsv_file, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            uniprot_id = line.split()[0]
            tsv_uniprot_ids.add(uniprot_id)

    # Step 2: Extract UniProt IDs from the FASTA file
    fasta_uniprot_ids = set()
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                uniprot_id = line.split('|')[0][1:]
                if "-" in uniprot_id:
                    uniprot_id = uniprot_id.split("-")[0]
                fasta_uniprot_ids.add(uniprot_id)

    # Step 3: Find UniProt IDs that are in the TSV but not in the FASTA file
    missing_uniprot_ids = tsv_uniprot_ids - fasta_uniprot_ids

    # Output the missing UniProt IDs
    if missing_uniprot_ids:
        print("The following UniProt IDs are in the TSV file but not in the FASTA file:")
        for uniprot_id in missing_uniprot_ids:
            print(uniprot_id)
    else:
        print("All UniProt IDs in the TSV file are present in the FASTA file.")

if __name__ == '__main__':
    input_tsv = './kinsnps/subkinsnps_uid_subs_split.txt'
    output_fasta = './kinsnps/subkinsnps.fasta'

    if os.path.exists(output_fasta):
        user_input = input(f"File {output_fasta} exists. Do you want to re-create it or append to it? (r/a): ").strip().lower()
        if user_input == 'r':
            mode = 'w'
        elif user_input == 'a':
            mode = 'a'
        else:
            print("Invalid option. Exiting.")
            exit()
    else:
        mode = 'w'

    read_tsv_and_generate_fasta(input_tsv, output_fasta, mode)
    report_missings(input_tsv, output_fasta)