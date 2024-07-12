import os
import csv

# Define the paths
allelic_variants_dir = './data/allelic_variants/'
omim_ids_csv_path = './data/omim_ids.csv'
kinase_list_csv_path = './data/kinase_list.csv'
output_file_path = './data/merged_allelic_variants.tsv'

# Define the output file headers
headers = ['OMIM_ID', 'Kinase_Description', 'Uniprot_ID', 'Gene', 'Number', 'Phenotype', 'Mutation', 'SNP', 'gnomAD_SNP', 'ClinVar']

# Read the OMIM IDs and gene information
omim_id_to_gene = {}
with open(omim_ids_csv_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        omim_id_to_gene[row['mimNumber']] = row['gene']

# Read the kinase list and uniprot ID information
gene_to_uniprot_id = {}
with open(kinase_list_csv_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_to_uniprot_id[row['gene']] = row['uniprot_id']

# Function to parse a single TSV file
def parse_tsv_file(tsv_path):
    with open(tsv_path, 'r') as file:
        lines = file.readlines()
        if len(lines) < 6:
            return []  # Invalid file format

        omim_id = lines[0].split("-")[1].strip()
        
        # Find the kinase_description
        kinase_desc = ""
        for i, line in enumerate(lines):
            if line.startswith("Allelic Variants ("):
                kinase_desc = lines[i-1].strip()
                break
        
        gene = omim_id_to_gene.get(omim_id, "Unknown")
        uniprot_id = gene_to_uniprot_id.get(gene, "Unknown")
        
        variants = []
        
        # Read the allelic variants section
        variant_lines = lines[i+3:]  # Allelic variants section starts 3 lines after "Allelic Variants ("
        for line in variant_lines:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 6:
                    number, phenotype, mutation, snp, gnomad_snp, clinvar = parts
                    variants.append([omim_id, kinase_desc, uniprot_id, gene, number, phenotype, mutation, snp, gnomad_snp, clinvar])
        
        return variants

# Collect all variants from all files
all_variants = []
for filename in os.listdir(allelic_variants_dir):
    if filename.endswith('.tsv'):
        tsv_path = os.path.join(allelic_variants_dir, filename)
        all_variants.extend(parse_tsv_file(tsv_path))

# Write the merged data to the output TSV file
with open(output_file_path, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t')
    writer.writerow(headers)
    writer.writerows(all_variants)

print(f"Merged TSV file created at: {output_file_path}")