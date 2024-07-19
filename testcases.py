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
    check_substitution(data, results, uniprot_id="Q96GX5", from_aa= "GLU", to_aa= "ASP", fullseq_pos= 167, alignment_pos= 128)
    check_substitution(data, results, uniprot_id="O14578", from_aa= "ASP", to_aa= "VAL", fullseq_pos= 230, alignment_pos= 128)
    check_substitution(data, results, uniprot_id="P22694", from_aa= "GLY", to_aa= "ARG", fullseq_pos= 235, alignment_pos= 187)
    check_substitution(data, results, uniprot_id="P22694", from_aa= "HIS", to_aa= "ASN", fullseq_pos= 88,  alignment_pos=41)
    check_substitution(data, results, uniprot_id="O14757", from_aa= "ARG", to_aa= "GLN", fullseq_pos= 379, alignment_pos= 379)
    check_substitution(data, results, uniprot_id="P57059", from_aa= "GLY", to_aa= "SER", fullseq_pos= 636, alignment_pos= 636)
    check_substitution(data, results, uniprot_id="Q13555", from_aa= "THR", to_aa= "MET", fullseq_pos= 240, alignment_pos= 209)
    check_substitution(data, results, uniprot_id="O14936", from_aa= "ARG", to_aa= "LEU", fullseq_pos= 28,  alignment_pos=17)
    check_substitution(data, results, uniprot_id="O14936", from_aa= "ASP", to_aa= "GLY", fullseq_pos= 710, alignment_pos= 710)
    check_substitution(data, results, uniprot_id="P48730", from_aa= "HIS", to_aa= "ARG", fullseq_pos= 46,  alignment_pos=38)
    check_substitution(data, results, uniprot_id="P21802", from_aa= "ALA", to_aa= "THR", fullseq_pos= 628, alignment_pos= 180)
    check_substitution(data, results, uniprot_id="P21802", from_aa= "GLU", to_aa= "ALA", fullseq_pos= 565, alignment_pos= 74)
    check_substitution(data, results, uniprot_id="P22607", from_aa= "ARG", to_aa= "HIS", fullseq_pos= 621, alignment_pos= 122)
    check_substitution(data, results, uniprot_id="P22607", from_aa= "ASN", to_aa= "SER", fullseq_pos= 540, alignment_pos= 59)
    check_substitution(data, results, uniprot_id="P22607", from_aa= "ASP", to_aa= "ASN", fullseq_pos= 628, alignment_pos= 129)
    check_substitution(data, results, uniprot_id="P22607", from_aa= "GLN", to_aa= "ARG", fullseq_pos= 485, alignment_pos= 14)

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