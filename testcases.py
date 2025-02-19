import json

def check_substitution(data, results, uniprot_id, from_aa, to_aa, fullseq_pos, alignment_pos):
    found = False
    for entry in data:
        if entry["uniprot_id"] == uniprot_id:
            for sub in entry["substitutions"]:
                if sub["from"] == from_aa and sub["to"] == to_aa:
                    if sub["full_sequence_pos"] == fullseq_pos and sub["alignment_pos"] == alignment_pos:
                        results["passed"].append(f"{uniprot_id}: Passed for {from_aa} to {to_aa} at fullseq_pos {fullseq_pos} and alignment_pos {alignment_pos}")
                    else:
                        if sub["full_sequence_pos"] != fullseq_pos:
                            results["failed"].append(f"{uniprot_id}: Expected full_sequence_pos {fullseq_pos}, got {sub['full_sequence_pos']}")
                        if sub["alignment_pos"] != alignment_pos:
                            results["failed"].append(f"{uniprot_id}: Expected alignment_pos {alignment_pos}, got {sub['alignment_pos']}")
                    found = True
                    break
            break
    if not found:
        results["failed"].append(f"{uniprot_id}: Entry with specified substitutions not found for {from_aa} to {to_aa}")

def test_output_json(output_file_path):
    results = {"passed": [], "failed": []}

    with open(output_file_path, 'r') as json_file:
        data = json.load(json_file)

    # Test cases
    check_substitution(data, results, uniprot_id="Q16644", from_aa= "L", to_aa= "P", fullseq_pos= 173, alignment_pos= 126)
    check_substitution(data, results, uniprot_id="Q8IW41", from_aa= "G", to_aa= "V", fullseq_pos= 107, alignment_pos= 79)
    check_substitution(data, results, uniprot_id="P68400", from_aa= "R", to_aa= "Q", fullseq_pos= 47, alignment_pos= 9)
    check_substitution(data, results, uniprot_id="P68400", from_aa= "D", to_aa= "H", fullseq_pos= 156, alignment_pos= 119)    
    check_substitution(data, results, uniprot_id="P68400", from_aa= "K", to_aa= "R", fullseq_pos= 198, alignment_pos= 160)    
    check_substitution(data, results, uniprot_id="P43405", from_aa= "S", to_aa= "Y", fullseq_pos= 550, alignment_pos= 171)
    check_substitution(data, results, uniprot_id="P10721", from_aa= "R", to_aa= "G", fullseq_pos= 796 , alignment_pos= 123)    
    check_substitution(data, results, uniprot_id="P10721", from_aa= "E", to_aa= "K", fullseq_pos= 839 , alignment_pos= 163)    
            
    return results

if __name__ == "__main__":
    output_file_path = 'kinsnps_allinfo.json'
    results = test_output_json(output_file_path)
    
    if results["passed"]:
        print("Passed tests:")
        for test in results["passed"]:
            print(f"  - {test}")
    else:
        print("No tests passed.")

    if results["failed"]:
        print("Failed tests:")
        for test in results["failed"]:
            print(f"  - {test}")
    else:
        print("No tests failed.")
