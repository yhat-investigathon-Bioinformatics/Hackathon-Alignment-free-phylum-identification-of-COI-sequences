import pandas as pd
import os
import re

# --- CONFIGURATION ---
INPUT_FILE = 'curated_dataset_CO1.csv'  # Output of curate_dataset.py
OUTPUT_DIR = 'phylums_fasta'            # Output directory for the FASTA files

def clean_name(name):
    """
    Cleans the phylum name so it is a valid file name.
    Replaces spaces and special characters with underscores.
    """
    if pd.isna(name) or name == "":
        return "Unknown"
    # Keep only letters, digits and underscores
    clean = re.sub(r'[^a-zA-Z0-9_]', '_', str(name))
    return clean

def write_fastas_by_phylum():
    print(f"Reading dataset from: {INPUT_FILE}...")
    try:
        df = pd.read_csv(INPUT_FILE)
    except FileNotFoundError:
        print(f"Error: '{INPUT_FILE}' not found. Run curate_dataset.py first.")
        return

    # Create output directory if it does not exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Directory '{OUTPUT_DIR}' created.")

    # List of unique phyla
    phylums = df['Phylum'].unique()
    print(f"Found {len(phylums)} distinct phyla.")

    print("-" * 40)
    print(f"{'PHYLUM':<30} | {'SEQUENCES':<10} | {'OUTPUT FILE'}")
    print("-" * 40)

    for phylum in phylums:
        # Filter rows for this phylum
        group = df[df['Phylum'] == phylum]
        count = len(group)

        # Skip empty groups (everything with > 0 sequences is written)
        if count == 0: continue

        # Build a safe file name
        safe_name = clean_name(phylum)
        filename = f"phylum_{safe_name}.fasta"
        filepath = os.path.join(OUTPUT_DIR, filename)

        # Write the FASTA file
        with open(filepath, 'w') as f:
            for index, row in group.iterrows():
                # FASTA format:
                # >SEQUENCE_ID
                # SEQUENCE
                f.write(f">{row['SequenceID']}\n")
                f.write(f"{row['Sequence']}\n")

        print(f"{str(phylum):<30} | {count:<10} | {filename}")

    print("-" * 40)
    print(f"Done. Output files are in '{OUTPUT_DIR}'.")

if __name__ == "__main__":
    write_fastas_by_phylum()
