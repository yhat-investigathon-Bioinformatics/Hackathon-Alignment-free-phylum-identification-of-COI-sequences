import numpy as np
from Bio import AlignIO
import math
import matplotlib.pyplot as plt
import sys

# --- CONFIGURATION ---
# Name of the aligned input file (MAFFT output)
ALIGNED_FILE = "aligned_Tubulinea.fasta"
# Length of the subsequence to find (e.g. 100 nucleotides)
WINDOW_SIZE = 100
OUTPUT_FILE = "optimal_subsequence.fasta"

def shannon_entropy(alignment):
    """Computes Shannon entropy (bits) for each alignment column."""
    print("Computing per-column entropy...")
    alignment_length = alignment.get_alignment_length()
    entropies = []

    # A numpy array would be faster for large datasets;
    # standard Biopython iteration is used for simplicity
    for i in range(alignment_length):
        column = alignment[:, i]
        col_length = len(column)

        # Count frequencies. Gaps could be excluded, but for separability
        # they may be informative (or noise).
        # Here the gap (-) is treated as just another character.
        counts = {}
        for base in column:
            if base not in counts: counts[base] = 0
            counts[base] += 1

        # Formula: H = -sum(p * log2(p))
        H = 0
        for base in counts:
            p = counts[base] / col_length
            if p > 0:
                H -= p * math.log2(p)

        entropies.append(H)
    return entropies

def sliding_window(entropy_values, window_size):
    """Finds the window with the highest mean entropy."""
    max_entropy = -1
    best_start = 0

    means = []

    # Iterate over the entropy values
    for i in range(len(entropy_values) - window_size + 1):
        # Mean entropy of the current window
        sub_window = entropy_values[i : i + window_size]
        current_mean = sum(sub_window) / window_size
        means.append(current_mean)

        if current_mean > max_entropy:
            max_entropy = current_mean
            best_start = i

    return best_start, max_entropy, means

def main():
    try:
        # 1. Load alignment
        print(f"Reading file: {ALIGNED_FILE}...")
        alignment = AlignIO.read(ALIGNED_FILE, "fasta")

        # 2. Compute entropy profile
        entropies = shannon_entropy(alignment)

        # 3. Find the best subsequence (sliding window)
        print(f"Searching for the best {WINDOW_SIZE} bp window...")
        best_start, max_ent, window_profile = sliding_window(entropies, WINDOW_SIZE)
        best_end = best_start + WINDOW_SIZE

        print("-" * 40)
        print("VARIABILITY ANALYSIS RESULT")
        print("-" * 40)
        print(f"Best subsequence found:")
        print(f" -> Start: position {best_start}")
        print(f" -> End:   position {best_end}")
        print(f" -> Mean entropy: {max_ent:.4f} bits")
        print("-" * 40)

        # 4. Save the subsequence to a new file
        # Trim the original alignment to keep only that region
        sub_alignment = alignment[:, best_start:best_end]
        AlignIO.write(sub_alignment, OUTPUT_FILE, "fasta")
        print(f"Subsequence extracted and saved to: {OUTPUT_FILE}")

        # 5. Plot for the report
        plt.figure(figsize=(12, 6))

        # Raw per-base entropy (light grey)
        plt.plot(entropies, color='lightgray', label='Per-base entropy')

        # Sliding-window mean (blue)
        # Shift the X axis so it aligns with the window center
        x_window = range(int(WINDOW_SIZE/2), int(WINDOW_SIZE/2) + len(window_profile))
        plt.plot(x_window, window_profile, color='blue', linewidth=2, label=f'Moving average ({WINDOW_SIZE}bp)')

        # Highlight the selected region
        plt.axvspan(best_start, best_end, color='red', alpha=0.3, label='Selected region')

        plt.title(f"Variability profile (Shannon) - {ALIGNED_FILE}")
        plt.xlabel("Alignment position")
        plt.ylabel("Entropy (bits)")
        plt.legend()
        plt.tight_layout()
        plt.show() # Or use plt.savefig("entropy_plot.png")

    except FileNotFoundError:
        print(f"ERROR: file '{ALIGNED_FILE}' not found.")
        print("Make sure the file name matches the MAFFT output.")

if __name__ == "__main__":
    main()
