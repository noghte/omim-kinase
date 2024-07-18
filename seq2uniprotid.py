import requests

def get_uniprot_id(protein_sequence):
    url = "https://www.uniprot.org/uniprot/"
    params = {
        'query': f'seq:{protein_sequence}',
        'format': 'tab',
        'columns': 'id'
    }
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        lines = response.text.splitlines()
        if len(lines) > 1:
            return lines[1].split()[0]  # The first line is the header, so we take the first entry of the second line
        else:
            return "No matching UniProt ID found."
    else:
        return f"Error: {response.status_code}"

# Example usage
protein_sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDPHSDEHGWQLVL"
uniprot_id = get_uniprot_id(protein_sequence)
print(f"UniProt ID: {uniprot_id}")