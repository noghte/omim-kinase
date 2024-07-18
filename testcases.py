import json

def check_substitution(data, uniprot_id, from_aa, to_aa, fullseq_pos, alignment_pos, results):
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
    check_substitution(data, "Q96GX5", "GLU", "ASP", 167, 128, results)
    check_substitution(data, "O14578", "ASP", "VAL", 230, 128, results)

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