import pandas as pd
import numpy as np
import re
import random
from Bio import SeqIO
from tqdm import tqdm
from collections import defaultdict

# Input files
INPUT_FILES = {
    "taxon": "MIDORI2_UNIQ_NUC_GB268_CO1.taxon",
    "fasta": "MIDORI2_UNIQ_NUC_GB268_CO1.fasta"
}

# Output prefix
OUTPUT_PREFIX = "dataset_10k"

# Criteria
MIN_PHYLUM_SIZE = 10000
SAMPLES_PER_CLASS = 10000
TRAIN_RATIO = 0.8
RANDOM_SEED = 42 # Fixed seed so the same dataset can be compared across models

# Data cleaning
def clean_midori_phylum(lineage_string):
    """
    Transforms: 'k__Eukaryota...;p__Arthropoda_555280;c__...'
    Into: 'Arthropoda'
    """
    if pd.isna(lineage_string): return None
    levels = str(lineage_string).split(';')
    for level in levels:
        level = level.strip()
        if level.startswith("p__"):
            # Remove 'p__' prefix and numeric suffixes
            name = re.sub(r'_\d+$', '', level.replace("p__", ""))
            # Additional specific cleanup
            if "class_" in name: name = name.replace("class_", "")
            return name
    return None

def main():
    print(f"--- STARTING ETL: CLEANING AND SPLIT ---")

    # ---------------------------------------------------------
    # STEP 1: READ AND CLEAN TAXONOMY
    # ---------------------------------------------------------
    print("\n Cleaning taxonomy and grouping...")

    phylum_to_ids = defaultdict(list)
    id_to_clean_label = {}

    try:
        df = pd.read_csv(INPUT_FILES["taxon"], sep='\t', header=None, dtype=str)
        ids = df.iloc[:, 0].str.strip().values
        lineages = df.iloc[:, 1].values

        for seq_id, lin in tqdm(zip(ids, lineages), total=len(ids)):
            phylum_clean = clean_midori_phylum(lin)

            if phylum_clean:
                phylum_to_ids[phylum_clean].append(seq_id)
                id_to_clean_label[seq_id] = phylum_clean # Store "Arthropoda"

    # Error handling
    except Exception as e:
        print(f"Error reading taxon file: {e}")
        return

    # ---------------------------------------------------------
    # STEP 2: SELECTION, BALANCING AND SPLIT
    # ---------------------------------------------------------
    print(f"\n Selecting phyla > {MIN_PHYLUM_SIZE} and splitting...")
    random.seed(RANDOM_SEED)

    train_ids = set()
    test_ids = set()
    final_phyla = []

    # Train/test split sizes
    n_train = int(SAMPLES_PER_CLASS * TRAIN_RATIO)
    n_test = SAMPLES_PER_CLASS - n_train

    sorted_phyla = sorted(phylum_to_ids.keys())

    print(f"   {'PHYLUM':<20} | {'TRAIN':<6} | {'TEST':<6}")
    print("-" * 45)

    for phylum in sorted_phyla:
        id_list = phylum_to_ids[phylum]
        total = len(id_list)

        if total >= MIN_PHYLUM_SIZE:
            # Random sampling of sequences
            selected = random.sample(id_list, SAMPLES_PER_CLASS)

            # Train/test split
            set_train = selected[:n_train]
            set_test = selected[n_train:]

            train_ids.update(set_train)
            test_ids.update(set_test)
            final_phyla.append(phylum)

            print(f"  {phylum:<20} | {len(set_train):<6} | {len(set_test):<6}")

    print("-" * 45)
    print(f"   Total Train IDs: {len(train_ids)}")
    print(f"   Total Test IDs:  {len(test_ids)}")

    # ---------------------------------------------------------
    # STEP 3: WRITE CLEAN FILES
    # ---------------------------------------------------------
    print("\n Writing output files...")

    # Output format: "ID \t Arthropoda"
    print("   -> Writing clean .taxon...")
    with open(f"{OUTPUT_PREFIX}_train.taxon", "w") as f_tr, open(f"{OUTPUT_PREFIX}_test.taxon", "w") as f_te:
        for seq_id, label in id_to_clean_label.items():
            # Only write IDs selected in the sampling step
            if seq_id in train_ids:
                f_tr.write(f"{seq_id}\t{label}\n")
            elif seq_id in test_ids:
                f_te.write(f"{seq_id}\t{label}\n")

    del id_to_clean_label

    print("   -> Writing .fasta...")
    count_train, count_test = 0, 0

    # Write .fasta files
    with open(f"{OUTPUT_PREFIX}_train.fasta", "w") as f_tr, open(f"{OUTPUT_PREFIX}_test.fasta", "w") as f_te:
        for record in tqdm(SeqIO.parse(INPUT_FILES["fasta"], "fasta")):
            seq_id = record.id.strip()

            if seq_id in train_ids:
                SeqIO.write(record, f_tr, "fasta")
                count_train += 1
            elif seq_id in test_ids:
                SeqIO.write(record, f_te, "fasta")
                count_test += 1

    print("\n ETL finished. Output files ready.")

if __name__ == "__main__":
    main()
