import pandas as pd

# Load the CSV files
csv_file_1 = './data/omim_ids_found_with_uniprot.csv'
csv_file_2 = './data/omim_ids_notfound_with_uniprot.csv'
df1 = pd.read_csv(csv_file_1)
df2 = pd.read_csv(csv_file_2)

# Combine both CSV files into a single DataFrame
df_combined = pd.concat([df1, df2])

# Load the TSV file
tsv_file = './data/merged_allelic_variants.tsv'
df_tsv = pd.read_csv(tsv_file, sep='\t')

# Merge the combined dataframe with the TSV dataframe on 'gene' and 'Gene' columns
merged_df = pd.merge(df_combined, df_tsv, left_on='gene', right_on='Gene', how='inner')

# Check if the uniprot_id from CSV matches the Uniprot_ID from TSV
merged_df['uniprot_id_match'] = merged_df['uniprot_id'] == merged_df['Uniprot_ID']

# Filter the rows where uniprot_id_match is False
mismatch_df = merged_df[~merged_df['uniprot_id_match']]

# Select only the required columns
result_df = mismatch_df[['gene', 'uniprot_id', 'Uniprot_ID']]
result_df.columns = ['gene', 'csv_uniprot_id', 'tsv_uniprot_id']

# Drop duplicate rows
result_df = result_df.drop_duplicates()

# Save the result to a new CSV file
output_file = './data/uniprot_id_mismatch.csv'
result_df.to_csv(output_file, index=False)

print(f"Mismatched UniProt IDs have been saved to {output_file}")