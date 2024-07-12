import os
import requests
import csv
from dotenv import load_dotenv
import time
import random

load_dotenv()
api_key = os.getenv('OMIM_API_KEY')

# Function to download and save allelic variants TSV file
def save_allelic_variants(mim_number):
    url = f"https://omim.org/allelicVariants/{mim_number}?format=tsv&apiKey={api_key}"
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate, br',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1'
    }
    
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        directory = './data/allelic_variants/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        filename = f"{directory}{mim_number}.tsv"
        with open(filename, 'w') as f:
            f.write(response.text)
        return filename
    elif response.status_code == 403:
        print(f"Access forbidden for MIM number {mim_number}: {response.text}")
    else:
        print(f"Failed to retrieve data for MIM number {mim_number}: Status code {response.status_code}")
    return None

if __name__ == "__main__":
    # Read the omim_ids.csv file
    with open('./data/omim_ids.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        omim_entries = [row for row in reader]

    # Save TSV files
    for row in omim_entries:
        if row["mimNumber"] != "Not Found":  # Ensure mimNumber is valid
            time.sleep(random.uniform(1, 3))  # Randomized delay between 1 and 3 seconds
            save_allelic_variants(row["mimNumber"])