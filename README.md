
## Preequistes

1. Register at [https://omim.org/api](https://omim.org/api) to get your API key
2. Rename the `.env.example` file to `.env` and set your API key in it

## Workflow to get omim data
**A)** Run `get_omim_ids.py` script to get omim ids

| It will iteate over the list of kinases that are in `data/kinase_list.csv` and get the omim ids. The output will be in `data/omim_ids.csv`

**B)** Run `save_allelic_variants_tsv_files.py` script to download allelic variants. 

|    It iterates over rows in `omim_ids.csv` and gets the allelic variants for each of them (e.g., if omim id is 131550, it will download `https://omim.org/allelicVariants/131550?format=tsv`). The output will be saved in `data/allelic_variants`

**C)** Run `get_uniprot_ids.py` script to add a column `uniprot_id` to the `omim_ids.csv`. 

**D)** Run `merge_tsvs.py` script to merge the allelic variants for each of the omim ids. The output will be in `data/merged_allelic_variants.tsv`   

**E)** With some Unix commands the file `subkinsnps_uid_subs_split.txt` has been created.

## Create JSON for analysis

Run `kinsnps_allinfo.py` script. Based on the files in `kinsnps` directory, it will produce `kinsnps_allinfo.json`.

---

## Other files

- `extract_sequences_from_uniprot.py`: Extract sequences of different isoforms from Uniprot.
  - Input:'./kinsnps/subkinsnps_uid_subs_split.txt'
  - Output: './kinsnps/subkinsnps.fasta' that needs to be aligned and the an `mma` file will be created.

- Then the `mma` file and the `subkinsnps_uid_subs_split.txt` should be used to create the `kinsnps_allinfo.json` file.
