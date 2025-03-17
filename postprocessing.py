import json
import os
import argparse
from typing import Dict, List

def load_json_file(file_path: str) -> List[Dict]:
    """Load data from JSON file."""
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: Could not find file {file_path}")
        return []
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in file {file_path}")
        return []

def save_json_file(data: List[Dict], file_path: str) -> None:
    """Save data to JSON file."""
    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)
        print(f"Successfully saved data to {file_path}")
    except Exception as e:
        print(f"Error saving data to {file_path}: {e}")

def cleanup_mutations(data: List[Dict]) -> List[Dict]:
    """Remove mutations where position exceeds sequence length."""
    cleaned_data = []
    removed_count = 0
    
    for protein in data:
        sequence = protein.get("sequence", "")
        seq_length = len([c for c in sequence if c not in '() -'])
        
        valid_substitutions = []
        for sub in protein.get("substitutions", []):
            pos = sub.get("full_sequence_pos", 0)
            if pos <= seq_length:
                valid_substitutions.append(sub)
            else:
                removed_count += 1
                
        # Create a new protein entry with valid substitutions only
        cleaned_protein = protein.copy()
        cleaned_protein["substitutions"] = valid_substitutions
        cleaned_data.append(cleaned_protein)
    
    print(f"Removed {removed_count} mutations with positions exceeding sequence length")
    return cleaned_data

def show_menu():
    """Display the main menu options."""
    print("\n=== Post-processing Menu ===")
    print("1. Cleanup mutations (remove positions exceeding sequence length)")
    print("0. Exit")
    return input("Select an option: ")

if __name__ == "__main__":
    """Main function to run the post-processing menu."""
    parser = argparse.ArgumentParser(description="Post-process kinase mutation data")
    parser.add_argument('--input', '-i', type=str, 
                      default="./data/kinsnps_allinfo_validmutations.json",
                      help='Path to input JSON file')
    parser.add_argument('--output', '-o', type=str, 
                      default=None,
                      help='Path to output JSON file (defaults to input file with _cleaned suffix)')
    
    args = parser.parse_args()
    input_file = args.input
    
    # Set default output file if not specified
    if args.output is None:
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_cleaned{ext}"
    else:
        output_file = args.output
    
    print(f"Using input file: {input_file}")
    print(f"Using output file: {output_file}")
    
    data = load_json_file(input_file)
    if not data:
        print("No data loaded. Exiting.")
        exit(1)
    
    while True:
        choice = show_menu()
        
        if choice == '1':
            print("\nCleaning up mutations...")
            cleaned_data = cleanup_mutations(data)
            save_json_file(cleaned_data, output_file)
        elif choice == '0':
            print("Exiting program.")
            break
        else:
            print("Invalid choice. Please try again.")