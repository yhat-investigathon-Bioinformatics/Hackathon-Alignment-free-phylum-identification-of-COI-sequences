import os
import sys
import subprocess
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO

# --- CONFIGURATION ---
# Input file
INPUT_FILE = "phylum_p__Discosea_555280.fasta"

# Analysis parameters
WINDOW_SIZE = 100  # Number of bases to extract
USE_FAST_MAFFT = False # Set to True for very large files (>5000 seqs)

def step_1_alignment(input_file):
    """
    Runs MAFFT from Python via a shell command.
    """
    base_name = os.path.splitext(input_file)[0]
    output_file = f"aligned_{base_name}.fasta"

    print(f"\n[STEP 1] Running MAFFT alignment for: {input_file}")
    print("         This may take a while depending on file size.")

    # Build the command
    if USE_FAST_MAFFT:
        # Fast mode for very large files
        cmd = f"mafft --retree 1 --maxiterate 0 {input_file} > {output_file}"
    else:
        # Automatic (accurate) mode
        cmd = f"mafft --auto {input_file} > {output_file}"

    # Run the command in the shell
    # shell=True allows using ">" to redirect output to a file
    process = subprocess.run(cmd, shell=True, executable='/bin/bash')

    if process.returncode != 0:
        print("Error: MAFFT failed. Make sure it is installed (sudo apt install mafft).")
        sys.exit(1)

    if not os.path.exists(output_file):
        print("Error: output file was not generated.")
        sys.exit(1)

    print(f"         Alignment completed: {output_file}")
    return output_file

def step_2_compute_entropy(aligned_file):
    """
    Computes Shannon entropy per column.
    """
    print(f"\n[STEP 2] Computing Shannon entropy (variability)...")
    try:
        alignment = AlignIO.read(aligned_file, "fasta")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        sys.exit(1)

    alignment_length = alignment.get_alignment_length()
    entropies = []

    # Iterate over columns
    for i in range(alignment_length):
        column = alignment[:, i]
        col_length = len(column)

        # Frequency count (gaps included as informative)
        counts = {}
        for base in column:
            counts[base] = counts.get(base, 0) + 1

        # Shannon: H = -sum(p * log2(p))
        H = 0
        for base in counts:
            p = counts[base] / col_length
            if p > 0:
                H -= p * math.log2(p)
        entropies.append(H)

    print(f"         Computation finished. Total length: {alignment_length} bp")
    return alignment, entropies

def step_3_extract_subsequence(alignment, entropies, source_file):
    """
    Finds the window with the highest mean entropy and saves it to a file.
    """
    print(f"\n[STEP 3] Searching for the most informative region (window: {WINDOW_SIZE} bp)...")

    max_mean = -1
    best_start = 0
    mean_profile = []

    # Sliding window
    # Note: this could be optimized with a convolution, but the loop is easier to read
    for i in range(len(entropies) - WINDOW_SIZE + 1):
        window = entropies[i : i + WINDOW_SIZE]
        mean = sum(window) / WINDOW_SIZE
        mean_profile.append(mean)

        if mean > max_mean:
            max_mean = mean
            best_start = i

    best_end = best_start + WINDOW_SIZE

    # Build the output file name
    base_name = os.path.splitext(source_file)[0]
    # Strip the "aligned_" prefix if present
    if base_name.startswith("aligned_"):
        base_name = base_name.replace("aligned_", "")

    output_sub = f"optimal_subsequence_{base_name}.fasta"

    # Trim and save
    sub_alignment = alignment[:, best_start:best_end]
    AlignIO.write(sub_alignment, output_sub, "fasta")

    print("-" * 50)
    print(f"FINAL RESULT FOR {source_file}")
    print("-" * 50)
    print(f"Best region found: position {best_start} to {best_end}")
    print(f"Mean entropy: {max_mean:.4f} bits")
    print(f"Output file: {output_sub}")
    print("-" * 50)

    return best_start, best_end, mean_profile, output_sub

def step_4_plot(entropies, mean_profile, start, end, file_name):
    """
    Generates the plot for the report.
    """
    print("\n[STEP 4] Generating plot...")
    plt.figure(figsize=(12, 6))

    # Raw data (grey)
    plt.plot(entropies, color='#cccccc', linewidth=0.5, label='Per-base entropy')

    # Moving average (blue)
    x_axis = range(int(WINDOW_SIZE/2), int(WINDOW_SIZE/2) + len(mean_profile))
    plt.plot(x_axis, mean_profile, color='blue', linewidth=2, label='Moving average')

    # Selected region (red)
    plt.axvspan(start, end, color='red', alpha=0.2, label='Optimal subsequence')

    plt.title(f"Variability analysis: {file_name}")
    plt.xlabel("Position")
    plt.ylabel("Entropy (bits)")
    plt.legend()

    # Save plot
    img_name = file_name.replace(".fasta", ".png")
    plt.savefig(img_name)
    print(f"         Plot saved as: {img_name}")
    # plt.show() # Uncomment if a graphical display is configured (e.g. in WSL)

# --- PIPELINE EXECUTION ---
if __name__ == "__main__":
    if not os.path.exists(INPUT_FILE):
        print(f"Error: file '{INPUT_FILE}' not found")
    else:
        # 1. Align
        aln_file = step_1_alignment(INPUT_FILE)

        # 2. Compute entropy
        aln_obj, entropy_list = step_2_compute_entropy(aln_file)

        # 3. Extract best region
        start, end, means, final_file = step_3_extract_subsequence(aln_obj, entropy_list, aln_file)

        # 4. Plot
        step_4_plot(entropy_list, means, start, end, final_file)

        print("\nPipeline finished successfully.")
