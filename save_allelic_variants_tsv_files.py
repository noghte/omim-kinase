import os
import requests
import csv
from dotenv import load_dotenv

load_dotenv()
api_key = os.getenv('OMIM_API_KEY')
# Function to download and save allelic variants TSV file
def save_allelic_variants(mim_number):
    url = f"https://omim.org/allelicVariants/{mim_number}?format=tsv&apiKey={api_key}"
    response = requests.get(url)
    if response.status_code == 200:
        directory = './data/allelic_variants/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        filename = f"{directory}{mim_number}.tsv"
        with open(filename, 'w') as f:
            f.write(response.text)
        return filename
    return None

if __name__ == "__main__":
    # Read the omim_ids.csv file
    with open('./data/omim_ids.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        omim_entries = [row for row in reader]

    # Save TSV files
    for row in omim_entries:
        if row["mimNumber"] != "Not Found":  # Ensure mimNumber is valid
            save_allelic_variants(row["mimNumber"])