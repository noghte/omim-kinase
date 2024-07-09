import os
import requests
import csv
from dotenv import load_dotenv

load_dotenv()
api_key = os.getenv('OMIM_API_KEY')

# Function to get OMIM ID and preferred title from the API
def get_omim_id_and_title(gene):
    url = f"https://api.omim.org/api/entry/search?search={gene}&exclude=all&format=json&start=0&limit=10&apiKey={api_key}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if data['omim']['searchResponse']['entryList']:
            entry = data['omim']['searchResponse']['entryList'][0]['entry']
            return entry['mimNumber'], entry['titles']['preferredTitle']
    return None, None

# Read the kinase list
kinase_list_file = './data/kinase_list.csv'
with open(kinase_list_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    kinases = [row['gene'] for row in reader]

# Write the OMIM IDs and titles to omim_ids.csv
output_file = './data/omim_ids.csv'
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['gene', 'mimNumber', 'preferredTitle']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for gene in kinases:
        mimNumber, preferredTitle = get_omim_id_and_title(gene)
        if mimNumber and preferredTitle:
            writer.writerow({'gene': gene, 'mimNumber': mimNumber, 'preferredTitle': preferredTitle})
        else:
            writer.writerow({'gene': gene, 'mimNumber': 'Not Found', 'preferredTitle': 'Not Found'})

print(f"OMIM IDs and titles have been written to {output_file}")