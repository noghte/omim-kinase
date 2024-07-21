import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
import random

# Load the CSV file
input_file = './data/omim_ids_notfound.csv'
df = pd.read_csv(input_file)

# Initialize the list to store UniProt IDs
uniprot_ids = []

# Output file path
output_file = './data/omim_ids_notfound_with_uniprot.csv'

# Create a copy of the dataframe to store the results
output_df = pd.DataFrame(columns=df.columns.tolist() + ['uniprot_id'])

# Iterate over the rows of the dataframe
for index, row in df.iterrows():
    mim_number = row['mimNumber']
    url = f'https://omim.org/entry/{mim_number}'
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate, br',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1'
    }
    try:
        # Fetch the HTML content
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Check if the request was successful
        html_content = response.text
        
        # Parse the HTML content
        soup = BeautifulSoup(html_content, 'html.parser')
        
        # Find the UniProt ID
        uniprot_link = soup.find('a', href=lambda href: href and 'uniprot.org/uniprotkb/' in href)
        if uniprot_link:
            uniprot_id = uniprot_link['href'].split('/')[-1]
        else:
            uniprot_id = None
        
    except requests.RequestException as e:
        print(f"Failed to fetch data for MIM number {mim_number}: {e}")
        uniprot_id = None
    
    # Append the UniProt ID to the list
    uniprot_ids.append(uniprot_id)
    
    # Create a new row with the UniProt ID
    new_row = row.to_dict()
    new_row['uniprot_id'] = uniprot_id
    
    # Append the new row to the output dataframe
    output_df = pd.concat([output_df, pd.DataFrame([new_row])], ignore_index=True)
    
    # Save the updated dataframe to the CSV file after each iteration
    output_df.to_csv(output_file, index=False)
    time.sleep(random.uniform(1, 4))  # Randomized delay between 1 and 4 seconds

print(f"Updated CSV file has been saved to {output_file}")