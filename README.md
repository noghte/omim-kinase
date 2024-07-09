
## Preequistes

1. Register at [https://omim.org/api](https://omim.org/api) to get your API key
2. Rename the `.env.example` file to `.env` and set your API key in it

## Workflow
**A)** Run `get_omim_ids.py` script to get omim ids

| It will iteate over the list of kinases that are in `kinase_list.csv` and get the omim ids. The output will be in `omim_ids.csv`

**B)** Run `get_allelic_variants.py` script to get allelic variants. 

|    It iterates over rows in `omim_ids.csv` and gets the allelic variants for each of them (e.g., if omim id is 131550, it will get data from `https://omim.org/allelicVariants/131550?format=tsv`). The output will be in `allelic_variants.csv`
   