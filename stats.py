import os
import csv

input_csv_path = './data/omim_ids.csv'
allelic_variants_dir = './data/allelic_variants/'
found_csv_path = './data/omim_ids_found.csv'
notfound_csv_path = './data/omim_ids_notfound.csv'

# Read the original CSV file
with open(input_csv_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    omim_entries = [row for row in reader]

# Separate entries based on the presence of the corresponding TSV file
found_entries = []
notfound_entries = []

for row in omim_entries:
    mim_number = row['mimNumber']
    tsv_filename = f"{allelic_variants_dir}{mim_number}.tsv"
    if os.path.exists(tsv_filename):
        found_entries.append(row)
    else:
        notfound_entries.append(row)

# Write found entries to omim_ids_found.csv
with open(found_csv_path, 'w', newline='') as csvfile:
    fieldnames = omim_entries[0].keys()
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(found_entries)

# Write not found entries to omim_ids_notfound.csv
with open(notfound_csv_path, 'w', newline='') as csvfile:
    fieldnames = omim_entries[0].keys()
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(notfound_entries)

print(f"Created {found_csv_path} with {len(found_entries)} entries.")
print(f"Created {notfound_csv_path} with {len(notfound_entries)} entries.")