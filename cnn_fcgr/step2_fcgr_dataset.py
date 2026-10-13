import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from tqdm import tqdm
import os

# ==========================================
# 1. GLOBAL CONFIGURATION
# ==========================================
# Input files
FILES = {
    "train_taxon": "dataset_10k_train.taxon",
    "train_fasta": "dataset_10k_train.fasta",
    "test_taxon":  "dataset_10k_test.taxon",
    "test_fasta":  "dataset_10k_test.fasta"
}

# Parameters (identical for Train and Test)
CENTER_POS = 225        # Detected entropy peak
WINDOW_LEN = 150        # Window size
SHIFTS = [-10, 0, 10]   # Data augmentation
K_MER = 6               # 64x64 pixels

# ==========================================
# 2. FCGR COMPUTATION
# ==========================================
def sequence_to_fcgr(sequence, k):
    """Converts a DNA sequence into a Frequency Chaos Game Representation image."""
    sequence = sequence.upper()
    table = str.maketrans("ACGT", "0123")
    s_clean = re.sub(r'[^ACGT]', '', sequence)

    if len(s_clean) < k: return None

    s_nums = [int(c) for c in s_clean.translate(table)]
    dim = 2**k
    matrix = np.zeros((dim, dim), dtype=np.float32)
    powers = [2**(k - 1 - j) for j in range(k)]

    for i in range(len(s_nums) - k + 1):
        kmer_nums = s_nums[i:i+k]
        row, col = 0, 0
        for j, base in enumerate(kmer_nums):
            p = powers[j]
            if base == 2: row += p; col += p # G
            elif base == 3: row += p         # T
            elif base == 1: col += p         # C
        matrix[row, col] += 1

    matrix = np.log1p(matrix)
    if np.max(matrix) > 0: matrix = matrix / np.max(matrix)
    return matrix

# ==========================================
# 3. PROCESSING
# ==========================================
def load_taxonomy(path):
    print(f"Loading taxonomy: {path}...")
    try:
        df = pd.read_csv(path, sep='\t', header=None, names=['ID', 'Phylum'])
        return pd.Series(df.Phylum.values, index=df.ID.str.strip()).to_dict()
    except Exception as e:
        print(f"Error loading taxonomy: {e}")
        return {}

def build_dataset(fasta_path, taxon_path, set_name):
    """Builds X and y arrays from the raw FASTA and taxonomy files."""
    print(f"\nPROCESSING {set_name.upper()} SET...")

    tax_map = load_taxonomy(taxon_path)
    if not tax_map: return None, None

    X_list = []
    y_list = []

    max_end_needed = CENTER_POS + max(SHIFTS) + WINDOW_LEN

    print(f"Reading sequences and generating FCGR images...")
    count = 0
    for record in tqdm(SeqIO.parse(fasta_path, "fasta")):
        seq_id = record.id.strip()
        phylum = tax_map.get(seq_id)

        # Skip records without a known label
        if not phylum: continue

        full_seq = str(record.seq).upper()
        if len(full_seq) < max_end_needed: continue

        for shift in SHIFTS:
            start = CENTER_POS + shift
            roi = full_seq[start : start + WINDOW_LEN]

            if roi.count('N') > (WINDOW_LEN * 0.1): continue

            img = sequence_to_fcgr(roi, K_MER)
            if img is not None:
                X_list.append(img)
                y_list.append(phylum)
                count += 1

    X = np.array(X_list, dtype=np.float32)
    y = np.array(y_list)

    if len(X) > 0:
        X = np.expand_dims(X, axis=-1)

    print(f"{set_name.upper()} generated: {len(X)} images.")
    return X, y

def plot_quality_control(X, y):
    """Generates the visual verification plot."""
    print("\nGenerating example control plot...")
    total_imgs = len(X)
    if total_imgs == 0: return

    indices = [0, total_imgs // 2, total_imgs - 1]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    cmap = 'inferno'

    for i, idx in enumerate(indices):
        ax = axes[i]
        img_2d = X[idx].squeeze()

        im = ax.imshow(img_2d, cmap=cmap, origin='upper')
        ax.set_title(f"Phylum: {y[idx]}", fontsize=12, fontweight='bold')
        ax.set_xlabel("K-mer coord. X", fontsize=9)
        if i == 0: ax.set_ylabel("K-mer coord. Y", fontsize=9)
        ax.set_xticks([0, 32, 63])
        ax.set_yticks([0, 32, 63])

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
    cbar.set_label('Normalized frequency (log)', fontsize=10)
    plt.suptitle(f"FCGR control: COI region centered at {CENTER_POS}bp", fontsize=14)

    output_img = "dataset_fcgr_control.png"
    plt.savefig(output_img, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {output_img}")
    plt.close()

# ==========================================
# 4. MAIN
# ==========================================
def main():
    print("=== FCGR DATASET GENERATOR (TRAIN & TEST) ===")

    # --- STEP 1: GENERATE TRAIN ---
    X_train, y_train = build_dataset(FILES["train_fasta"], FILES["train_taxon"], "Train")

    if X_train is not None and len(X_train) > 0:
        # Save Train
        np.save("X_train.npy", X_train)
        np.save("y_train.npy", y_train)
        print("X_train.npy and y_train.npy saved.")

        # Visualize (Train only, for examples)
        plot_quality_control(X_train, y_train)
    else:
        print("No Train data generated. Check the input files.")

    # --- STEP 2: GENERATE TEST ---
    X_test, y_test = build_dataset(FILES["test_fasta"], FILES["test_taxon"], "Test")

    if X_test is not None and len(X_test) > 0:
        # Save Test
        np.save("X_test.npy", X_test)
        np.save("y_test.npy", y_test)
        print("X_test.npy and y_test.npy saved.")
    else:
        print("No Test data generated. Check the input files.")

    print("\n=== PROCESS COMPLETED ===")

if __name__ == "__main__":
    main()
