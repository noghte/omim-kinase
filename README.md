## Preequistes

1. Register at [https://omim.org/api](https://omim.org/api) to get your API key
2. Rename the `.env.example` file to `.env` and set your API key in it

## Workflow to get OMIM data

### Prerequisites
1. Register at [https://omim.org/api](https://omim.org/api) to get your API key
2. Rename the `.env.example` file to `.env` and set your API key in it

### Step 1: Get OMIM IDs
Run `get_omim_ids.py` to fetch OMIM IDs for your kinases:
- **Input**: `data/kinase_list.csv` (list of kinase genes)
- **Output**: `data/omim_ids.csv` with columns:
  - `gene`: Kinase gene name
  - `mimNumber`: OMIM ID
  - `preferredTitle`: OMIM entry title
- **Reference**: 

### Step 2: Download Allelic Variants
Run `save_allelic_variants_tsv_files.py` to fetch variant data:
- **Input**: `data/omim_ids.csv`
- **Output**: TSV files in `data/allelic_variants/` directory (one file per OMIM ID)
- **Format**: Each TSV contains:
  - Header information
  - Variant table with columns: Number, Phenotype, Mutation, SNP, gnomAD_SNP, ClinVar
- **Note**: Includes rate limiting (1-3 second delay between requests)
- **Reference**:

### Step 3: Add UniProt IDs
Run `get_uniprot_ids.py` to add UniProt mappings:
- **Input**: `data/omim_ids.csv`
- **Output**: Updated CSV with new `uniprot_id` column
- **Process**: Scrapes UniProt IDs from OMIM entry pages
- **Reference**:

### Step 4: Merge Variant Data
Run `merge_tsvs.py` to combine all variant information:
- **Input**: 
  - `data/allelic_variants/*.tsv` files
  - `data/omim_ids.csv`
  - `data/kinase_list.csv`
- **Output**: `data/merged_allelic_variants.tsv`
- **Columns**: OMIM_ID, Kinase_Description, Uniprot_ID, Gene, Number, Phenotype, Mutation, SNP, gnomAD_SNP, ClinVar
- **Reference**:

### Step 5: Generate Final Output
Use Unix commands to process `data/merged_allelic_variants.tsv` and create `subkinsnps_uid_subs_split.txt`

### Optional: Generate Statistics
Run `stats.py` to check data completeness:
- **Output**: 
  - `data/omim_ids_found.csv`: Successfully processed entries
  - `data/omim_ids_notfound.csv`: Failed/missing entries
- **Reference**:

## Workflow to get ClinVar data

**A)** Run `get_clinvar.py`. It iterates over `omim_ids.csv` and calls https://www.ncbi.nlm.nih.gov/clinvar?term={mimNumber}[MIM] to get the ClinVar data. It clicks on the "Create File" button and downloads the text file and saves it to `./data/clinvar/{mimNumber}.txt`.


## Create JSON for analysis

Run `kinsnps_allinfo.py` script. Based on the files in `kinsnps` directory, it will produce `kinsnps_allinfo.json`.

---

## Other files

- `extract_sequences_from_uniprot.py`: Extract sequences of different isoforms from Uniprot.
  - Input:'./kinsnps/subkinsnps_uid_subs_split.txt'
  - Output: './kinsnps/subkinsnps.fasta' that needs to be aligned and the an `mma` file will be created.

- Then the `mma` file and the `subkinsnps_uid_subs_split.txt` should be used to create the `kinsnps_allinfo.json` file.
