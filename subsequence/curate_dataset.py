import pandas as pd
from Bio import SeqIO
import sys

# --- FILE CONFIGURATION ---
# Adjust these names if your files use a different extension (e.g. .txt or .tsv)
FILE_FASTA = 'MIDORI2_UNIQ_NUC_GB268_CO1.fasta'
FILE_TAXON = 'MIDORI2_UNIQ_NUC_GB268_CO1.taxon'
FILE_OUTPUT = 'curated_dataset_CO1.csv'

def load_taxonomy(taxon_path):
    """
    Loads the taxonomy file into a dictionary.
    Maps SEQUENCE_ID -> TAXONOMY_STRING
    """
    print(f"Loading taxonomy from {taxon_path}...")
    tax_map = {}
    try:
        # Plain open() instead of gzip.open()
        with open(taxon_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    # parts[0] is the ID, parts[1] is the full taxonomy
                    tax_map[parts[0]] = parts[1]
    except FileNotFoundError:
        print(f"Error: file '{taxon_path}' not found. Check the file name.")
        sys.exit(1)
    except UnicodeDecodeError:
        print("Encoding error: try changing encoding='utf-8' to 'latin-1'.")
        sys.exit(1)

    print(f"-> Taxonomy loaded: {len(tax_map)} entries.")
    return tax_map

def extract_phylum(tax_string):
    """
    Extracts the phylum from the taxonomy string.
    Typical structure: Eukaryota;Arthropoda;Insecta...
    """
    if not tax_string: return "Unknown"
    parts = tax_string.split(';')
    # Assume the phylum is the second element (index 1)
    if len(parts) > 1:
        return parts[1]
    return parts[0]

def process_data():
    # 1. Load the taxonomy map first
    tax_map = load_taxonomy(FILE_TAXON)

    data_list = []
    valid_bases = set("ACGT") # Strict set used to validate sequences

    print(f"Processing FASTA from {FILE_FASTA}...")

    count_total = 0
    count_length_discard = 0
    count_ambiguous_discard = 0
    count_kept = 0

    try:
        # Plain file open
        with open(FILE_FASTA, "r", encoding='utf-8') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                count_total += 1

                # Uppercase sequence
                seq_str = str(record.seq).upper()
                length = len(seq_str)

                # --- FILTER 1: LENGTH [300, 1700] ---
                if not (300 <= length <= 1700):
                    count_length_discard += 1
                    continue

                # --- FILTER 2: AMBIGUOUS BASES ---
                # Keep only sequences containing A, C, G, T
                if not set(seq_str).issubset(valid_bases):
                    count_ambiguous_discard += 1
                    continue

                # --- TAXONOMY LOOKUP ---
                # Clean the ID (keep the part before the first space)
                seq_id = record.id.split()[0]

                # Look up the taxonomy
                full_tax = tax_map.get(seq_id, "Unknown")

                # Extract the phylum
                phylum = extract_phylum(full_tax)

                # Store the record
                data_list.append({
                    'SequenceID': seq_id,
                    'Phylum': phylum,
                    'Taxonomy': full_tax,
                    'Sequence': seq_str,
                    'Length': length
                })

                count_kept += 1

                # Progress report every 50,000 sequences
                if count_total % 50000 == 0:
                    print(f" Processed: {count_total} | Kept: {count_kept}...")

    except FileNotFoundError:
        print(f"Error: file '{FILE_FASTA}' not found.")
        sys.exit(1)

    # --- SAVE RESULTS ---
    print("\nBuilding DataFrame and saving file...")
    df = pd.DataFrame(data_list)

    # Save as plain CSV (uncompressed)
    df.to_csv(FILE_OUTPUT, index=False)

    print("-" * 30)
    print("FINAL REPORT:")
    print(f"Total read: {count_total}")
    print(f"Discarded by length (<300 or >1700): {count_length_discard}")
    print(f"Discarded by ambiguous bases (N, M, R...): {count_ambiguous_discard}")
    print(f"CURATED SEQUENCES KEPT: {count_kept}")
    print(f"File saved to: {FILE_OUTPUT}")
    print("-" * 30)

    return df

if __name__ == "__main__":
    df_result = process_data()
    print(df_result.head())
