from typing import Dict
import csv
import requests
from bs4 import BeautifulSoup
from time import sleep
import random

def get_all_omim_ids() -> set:
    """Extract all OMIM IDs from the TSV files."""
    omim_ids = set()
    # read from omim_ids.csv
    with open('./data/omim_ids.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['mimNumber'].isdigit():
                omim_ids.add(row['mimNumber'])
    return omim_ids

def create_omim_uniprot_mapping() -> Dict[str, str]:
    """Create mapping between OMIM IDs and Uniprot IDs."""
    mapping = {}
    
    # Process found mappings
    with open('./data/omim_uniprot_mapping.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            mapping[row['mimNumber']] = row['uniprot_id']
    
    return mapping

def get_uniprot_ids_from_fasta(fasta_file: str) -> set:
    """Extract UniProt IDs from FASTA file headers."""
    uniprot_ids = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract UniProt ID (everything before the first |)
                uniprot_id = line.strip().split('|')[0][1:]  # Remove '>' and get ID
                uniprot_ids.add(uniprot_id)
    return uniprot_ids

if __name__ == '__main__':
    omim_ids = get_all_omim_ids()
    omim_uniprot_mapping = {}
    for omim_id in list(omim_ids):
            url = f'https://omim.org/entry/{omim_id}'
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.5',
                'Accept-Encoding': 'gzip, deflate, br',
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1'
            }
            try:
                response = requests.get(url, headers=headers)
                sleep(random.uniform(0.5, 3.5))
                response.raise_for_status()
                html_content = response.text
                soup = BeautifulSoup(html_content, 'html.parser')
                uniprot_link = soup.find('a', href=lambda href: href and 'uniprot.org/uniprotkb/' in href)
                uniprot_id = uniprot_link['href'].split('/')[-1]
                omim_uniprot_mapping[omim_id] = uniprot_id
                print(f"OMIM ID: {omim_id}, UniProt ID: {uniprot_id}")
            except requests.RequestException as e:
                print(f"Failed to fetch data for omim id {omim_id}: {e}")
                uniprot_id = None
    #save
    with open('./data/omim_uniprot_mapping.csv', 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['omim_id', 'uniprot_id'])
        for omim_id, uniprot_id in omim_uniprot_mapping.items():
            writer.writerow([omim_id, uniprot_id])
    
    # check fasta file
    uniprot_data = create_omim_uniprot_mapping()

    input_fasta = './kinsnps/human_kinases.fasta'

    # Get UniProt IDs from FASTA file
    fasta_uniprot_ids = get_uniprot_ids_from_fasta(input_fasta)
    
    # Find missing IDs
    missing_ids = set(uniprot_data.values()) - fasta_uniprot_ids
    
    # Print results
    print(f"Total UniProt IDs in mapping: {len(uniprot_data)}")
    print(f"Total UniProt IDs in FASTA: {len(fasta_uniprot_ids)}")
    print(f"Number of missing IDs: {len(missing_ids)}")
    print("\nMissing UniProt IDs:")
    for missing_id in sorted(missing_ids):
        print(missing_id)
