from typing import Dict
import csv
import requests
from bs4 import BeautifulSoup
from time import sleep
import random
import os
from collections import OrderedDict

user_agents = [
    # Chrome (Windows/Mac/Linux)
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36',
    
    # Firefox
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:109.0) Gecko/20100101 Firefox/115.0',
    
    # Safari
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.1 Safari/605.1.15',
    
    # Edge
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36 Edg/121.0.0.0',
    
    # Mobile (iOS/Android)
    'Mozilla/5.0 (iPhone; CPU iPhone OS 17_1_1 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.1 Mobile/15E148 Safari/604.1'
]
def get_headers():
    return OrderedDict([
        ('User-Agent', random.choice(user_agents)),
        ('Accept', 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8'),
        ('Accept-Language', random.choice(['en-US,en;q=0.5', 'en-GB,en;q=0.7', 'en-CA,en;q=0.3'])),
        ('Accept-Encoding', 'gzip, deflate, br'),
        # ('Referer', random.choice(referers)),
        ('Connection', 'keep-alive'),
        ('Upgrade-Insecure-Requests', '1'),
        ('Sec-Fetch-Dest', 'document'),
        ('Sec-Fetch-Mode', 'navigate'),
        ('Sec-Fetch-Site', 'same-origin'),
        ('Sec-Fetch-User', '?1'),
        ('DNT', random.choice(['1', '0']))
    ])
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

def get_notmapped_omim_ids() -> set:
    """Get OMIM IDs that don't have UniProt mappings in the CSV file."""
    mapped_ids = set()
    all_omim_ids = get_all_omim_ids()
    
    # Read mapped IDs from CSV if it exists
    csv_file_path = './data/omim_uniprot_mapping.csv'
    if os.path.exists(csv_file_path):
        with open(csv_file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                mapped_ids.add(row['omim_id'])
    
    # Return the difference between all OMIM IDs and mapped IDs
    return all_omim_ids - mapped_ids
    
if __name__ == '__main__':
    # Check for unmapped OMIM IDs first
    unmapped_ids = get_notmapped_omim_ids()
    print(f"Found {len(unmapped_ids)} unmapped OMIM IDs:")
    for omim_id in sorted(unmapped_ids):
        print(omim_id)
    print("-" * 50)

    omim_ids = get_all_omim_ids()
    omim_uniprot_mapping = {}
    # 119530
    # 613762
    # 616731
    # 617104
    # Before the loop, check if CSV file exists and load existing mappings
    existing_mappings = {}
    csv_file_path = './data/omim_uniprot_mapping.csv'

    # Try to read existing mappings if the file exists
    if os.path.exists(csv_file_path):
        with open(csv_file_path, 'r', newline='') as csvfile:
            csv_reader = csv.reader(csvfile)
            next(csv_reader)  # Skip header row
            for row in csv_reader:
                if len(row) >= 2:
                    omim_id, uniprot_id = row[0], row[1]
                    existing_mappings[omim_id] = uniprot_id
                    omim_uniprot_mapping[omim_id] = uniprot_id  # Also update the main mapping
        print(f"Loaded {len(existing_mappings)} existing mappings from CSV file")

    # Now open the file in append mode if it exists, or write mode if it doesn't
    file_mode = 'a' if os.path.exists(csv_file_path) else 'w'
    with open(csv_file_path, file_mode, newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write header row only if we're creating a new file
        if file_mode == 'w':
            csv_writer.writerow(['OMIM_ID', 'UniProt_ID'])
        
        for omim_id in list(omim_ids):
            # Skip if we already have this mapping
            if omim_id in existing_mappings:
                print(f"Skipping OMIM ID: {omim_id}, already mapped to UniProt ID: {existing_mappings[omim_id]}")
                continue
                
            url = f'https://omim.org/entry/{omim_id}'
            headers = get_headers()
            try:
                response = requests.get(url, headers=headers)
                sleep(random.uniform(0.5, 2.5))
                response.raise_for_status()
                html_content = response.text
                soup = BeautifulSoup(html_content, 'html.parser')
                uniprot_link = soup.find('a', href=lambda href: href and 'uniprot.org/uniprotkb/' in href)
                uniprot_id = uniprot_link['href'].split('/')[-1]
                omim_uniprot_mapping[omim_id] = uniprot_id
                print(f"OMIM ID: {omim_id}, UniProt ID: {uniprot_id}")
                
                # Write to CSV
                csv_writer.writerow([omim_id, uniprot_id])
                
            except requests.RequestException as e:
                print(f"Failed to fetch data for omim id {omim_id}: {e}")
                uniprot_id = None
                
    print("Done!")


    # check fasta file
    # uniprot_data = create_omim_uniprot_mapping()

    # input_fasta = './kinsnps/human_kinases.fasta'

    # # Get UniProt IDs from FASTA file
    # fasta_uniprot_ids = get_uniprot_ids_from_fasta(input_fasta)
    
    # # Find missing IDs
    # missing_ids = set(uniprot_data.values()) - fasta_uniprot_ids
    
    # # Print results
    # print(f"Total UniProt IDs in mapping: {len(uniprot_data)}")
    # print(f"Total UniProt IDs in FASTA: {len(fasta_uniprot_ids)}")
    # print(f"Number of missing IDs: {len(missing_ids)}")
    # print("\nMissing UniProt IDs:")
    # for missing_id in sorted(missing_ids):
    #     print(missing_id)
