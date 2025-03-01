import json

def check_substitution(data, results, uniprot_id, from_aa, to_aa, fullseq_pos, alignment_pos):
    found = False
    for entry in data:
        if entry["uniprot_id"] == uniprot_id:
            for sub in entry["substitutions"]:
                if sub["from"] == from_aa and sub["to"] == to_aa:
                    if sub["full_sequence_pos"] == fullseq_pos and sub["alignment_pos"] == alignment_pos:
                        results["passed"].append(f"{uniprot_id}: Passed for {from_aa} to {to_aa} at fullseq_pos {fullseq_pos} and alignment_pos {alignment_pos}")
                        found = True
                        return  # Found exact match
                    elif not found:  # Only record failure if no match found yet
                        if sub["full_sequence_pos"] != fullseq_pos:
                            results["failed"].append(f"{uniprot_id}: Expected full_sequence_pos {fullseq_pos}, got {sub['full_sequence_pos']}")
                        if sub["alignment_pos"] != alignment_pos:
                            results["failed"].append(f"{uniprot_id}: Expected alignment_pos {alignment_pos}, got {sub['alignment_pos']}")
                        found = True  # Mark as found even if positions don't match
    
    if not found:
        results["failed"].append(f"{uniprot_id}: Entry with specified substitutions not found for {from_aa} to {to_aa}")

def test_output_json(json_file_path):
    results = {"passed": [], "failed": []}
    
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        results["failed"].append(f"Error reading JSON file: {str(e)}")
        return results

    check_substitution(data, results, "Q16644", "L", "P", 173, 126)
    check_substitution(data, results, "Q8IW41", "G", "V", 107, 79)
    check_substitution(data, results, "P68400", "R", "Q", 47, 9)
    check_substitution(data, results, "P68400", "D", "H", 156, 119)
    check_substitution(data, results, "P68400", "K", "R", 198, 160)
    check_substitution(data, results, "P43405", "S", "Y", 550, 171)
    check_substitution(data, results, "P10721", "R", "G", 796, 123)
    check_substitution(data, results, "P10721", "E", "K", 839, 163)

    return results

if __name__ == "__main__":
    results = test_output_json('kinsnps_allinfo.json')
    
    print("\nTest Results:")
    if results["passed"]:
        print("\nPassed tests:")
        for test in results["passed"]:
            print(f"  ✓ {test}")
            
    if results["failed"]:
        print("\nFailed tests:")
        for test in results["failed"]:
            print(f"  ✗ {test}")
            
    if not results["passed"] and not results["failed"]:
        print("No tests were executed")