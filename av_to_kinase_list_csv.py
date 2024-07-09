import json

rows = []
with open('./data/av.json', 'r') as f:
    data = json.load(f)["content"]
    for key, value in data.items():
        row = (value["UniProt_ID"],key)
        rows.append(row)

with open('./data/kinase_list.csv', 'w') as f:
    f.write('uniprot_id,gene\n')
    for row in rows:
        f.write(f'{row[0]},{row[1]}\n')
