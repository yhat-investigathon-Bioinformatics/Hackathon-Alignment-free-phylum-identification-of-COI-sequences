import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
from tqdm import tqdm
import math
from scipy.ndimage import gaussian_filter1d

# ==========================================
# 1. CONFIGURATION
# ==========================================
FILES = {
    "taxon": "dataset_10k_train.taxon",
    "fasta": "dataset_10k_train.fasta"
}

# The 8 largest phyla
TARGET_PHYLA = [
    "Annelida", "Arthropoda", "Chordata", "Echinodermata",
    "Mollusca", "Nematoda", "Platyhelminthes", "Rhodophyta"
]

# Plot colors
COLORS = {
    'Annelida': 'brown',
    'Arthropoda': 'grey',
    'Chordata': 'black',
    'Echinodermata': 'green',
    'Mollusca': 'red',
    'Nematoda': 'purple',
    'Platyhelminthes': 'blue',
    'Rhodophyta': 'magenta'
}

# Scan parameters
K_MER = 6
WINDOW_SIZE = 150
STEP_SIZE = 25
MAX_LEN_ANALYSIS = 1000

# ==========================================
# 2. FUNCTIONS
# ==========================================
def load_taxonomy(path):
    print(f"Loading taxonomy from {path}...")
    try:
        df = pd.read_csv(path, sep='\t', header=None, names=['ID', 'Phylum'])
        return pd.Series(df.Phylum.values, index=df.ID.str.strip()).to_dict()
    # Error handling
    except Exception as e:
        print(f"Error: {e}")
        return {}

# Shannon entropy
def shannon_entropy(kmers):
    total = len(kmers)
    if total == 0: return 0.0
    counts = Counter(kmers)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy

# ==========================================
# 3. ANALYSIS
# ==========================================
def main():
    tax_map = load_taxonomy(FILES["taxon"])
    if not tax_map: return

    profiles = {p: {} for p in TARGET_PHYLA}

    print(f"Computing entropy for the 8 phyla...")
    count = 0

    for record in tqdm(SeqIO.parse(FILES["fasta"], "fasta")):
        seq_id = record.id.strip()
        phylum = tax_map.get(seq_id)

        if phylum not in TARGET_PHYLA: continue

        seq = str(record.seq).upper()
        if len(seq) < 600: continue

        count += 1
        limit = min(len(seq), MAX_LEN_ANALYSIS)

        for start in range(0, limit - WINDOW_SIZE, STEP_SIZE):
            window = seq[start : start + WINDOW_SIZE]
            kmers = [window[i:i+K_MER] for i in range(len(window) - K_MER + 1)]
            h = shannon_entropy(kmers)

            if start not in profiles[phylum]:
                profiles[phylum][start] = []
            profiles[phylum][start].append(h)

    print(f"Processed {count} sequences.")

    # ==========================================
    # 4. PLOT
    # ==========================================
    TARGET = 225
    WIDTH = WINDOW_SIZE # 150

    # Compute bounds to center the window
    win_start = TARGET - (WIDTH // 2) # 225 - 75 = 150
    win_end = TARGET + (WIDTH // 2)    # 225 + 75 = 300

    print("Generating plot...")
    plt.figure(figsize=(14, 8))

    for phylum in TARGET_PHYLA:
        data = profiles[phylum]
        if not data: continue

        pos = sorted(data.keys())
        vals = [np.mean(data[p]) for p in pos]
        y_smooth = gaussian_filter1d(vals, sigma=1.5)

        lw = 3 if phylum in ['Mollusca', 'Platyhelminthes'] else 1.5
        alpha = 1.0 if phylum in ['Mollusca', 'Platyhelminthes'] else 0.6
        ls = '--' if phylum == 'Arthropoda' else '-'

        plt.plot(pos, y_smooth, label=phylum, color=COLORS.get(phylum, 'black'),
                 linewidth=lw, alpha=alpha, linestyle=ls)

    # Plot settings
    plt.axvline(x=TARGET, color='black', linestyle=':', label=f'Optimal center ({TARGET})')
    plt.axvspan(win_start, win_end, color='yellow', alpha=0.3, label=f'Window ({win_start}-{win_end})')

    plt.title(f"Shannon entropy by position (K={K_MER})", fontsize=16)
    plt.xlabel("Position (bp)", fontsize=12)
    plt.ylabel("Shannon entropy", fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("subseq_co1.png", dpi=150)
    print("Plot saved: subseq_co1.png")
    plt.show()

if __name__ == "__main__":
    main()
