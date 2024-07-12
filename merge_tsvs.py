import os
import csv

# Define the paths
allelic_variants_dir = './data/allelic_variants/'
output_file_path = './data/merged_allelic_variants.tsv'

# Define the output file headers
headers = ['OMIM_ID', 'Kinase_Name', 'Number', 'Phenotype', 'Mutation', 'SNP', 'gnomAD_SNP', 'ClinVar']

# Function to parse a single TSV file
def parse_tsv_file(tsv_path):
    with open(tsv_path, 'r') as file:
        lines = file.readlines()
        if len(lines) < 6:
            return []  # Invalid file format

        omim_id = lines[0].split("-")[1].strip()
        
        # Find the kinase_name
        kinase_name = ""
        for i, line in enumerate(lines):
            if line.startswith("Allelic Variants ("):
                kinase_name = lines[i-1].strip()
                break
        
        variants = []
        
        # Read the allelic variants section
        variant_lines = lines[i+3:]  # Allelic variants section starts 3 lines after "Allelic Variants ("
        for line in variant_lines:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 6:
                    number, phenotype, mutation, snp, gnomad_snp, clinvar = parts
                    variants.append([omim_id, kinase_name, number, phenotype, mutation, snp, gnomad_snp, clinvar])
        
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