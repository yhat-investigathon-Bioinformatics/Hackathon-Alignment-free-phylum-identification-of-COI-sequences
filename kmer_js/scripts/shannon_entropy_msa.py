#!/usr/bin/env python3
# Usage:
#   python3 shannon_entropy_msa.py pilot_300_all.aln.fasta 200 entropy_pilot
#
# Where:
#   - 200 = window size (in alignment columns)
#   - entropy_pilot = output prefix

import sys
import math

def read_alignment(fasta_path):
    """
    Read an aligned FASTA file and return:
      - list of sequences (strings)
      - alignment length (all sequences must have the same length)
    """
    seqs = []
    curr = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if curr:
                    seqs.append("".join(curr))
                    curr = []
            else:
                curr.append(line)
        if curr:
            seqs.append("".join(curr))

    if not seqs:
        raise ValueError("No sequences were read from the alignment")

    aln_len = len(seqs[0])
    for s in seqs:
        if len(s) != aln_len:
            raise ValueError("Sequences of different lengths in the MSA")

    return seqs, aln_len

def shannon_entropy_column(chars):
    """
    Compute Shannon entropy (in bits) for a column.
    Ignores gaps ('-' and '.').
    """
    counts = {}
    total = 0
    for c in chars:
        c = c.upper()
        if c in ("-", ".", "N"):  # treat gaps and N as non-informative
            continue
        counts[c] = counts.get(c, 0) + 1
        total += 1

    if total == 0 or len(counts) <= 1:
        return 0.0

    h = 0.0
    for v in counts.values():
        p = v / total
        h -= p * math.log2(p)
    return h

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 shannon_entropy_msa.py aln.fasta window_size out_prefix")
        sys.exit(1)

    aln_file = sys.argv[1]
    window_size = int(sys.argv[2])
    out_prefix = sys.argv[3]

    print(f"[INFO] Reading alignment from {aln_file} ...", file=sys.stderr)
    seqs, aln_len = read_alignment(aln_file)
    nseq = len(seqs)
    print(f"[INFO] nseq = {nseq}, aln_len = {aln_len}", file=sys.stderr)

    # --- Per-column entropy ---
    entropies = []
    gap_fracs = []

    print("[INFO] Computing per-column entropy...", file=sys.stderr)
    for pos in range(aln_len):
        col = [s[pos] for s in seqs]
        gap_count = sum(1 for c in col if c in ("-", "."))
        gap_frac = gap_count / nseq
        H = shannon_entropy_column(col)
        entropies.append(H)
        gap_fracs.append(gap_frac)

    # Save per-column TSV
    col_tsv = out_prefix + "_per_column.tsv"
    with open(col_tsv, "w") as out:
        out.write("col_idx\tentropy_bits\tgap_fraction\n")
        for i, (H, g) in enumerate(zip(entropies, gap_fracs), start=1):
            out.write(f"{i}\t{H:.6f}\t{g:.6f}\n")

    print(f"[INFO] Per-column entropy saved to {col_tsv}", file=sys.stderr)

    # --- Sliding-window entropy ---
    if window_size > aln_len:
        print("[WARN] window_size > aln_len; skipping window computation", file=sys.stderr)
        sys.exit(0)

    win_tsv = out_prefix + f"_window{window_size}.tsv"
    print(f"[INFO] Computing mean entropy over windows of size {window_size}...", file=sys.stderr)

    with open(win_tsv, "w") as out:
        out.write("window_start\twindow_end\tmean_entropy_bits\tmean_gap_fraction\n")
        for start in range(0, aln_len - window_size + 1):
            end = start + window_size  # end is exclusive
            window_H = entropies[start:end]
            window_gap = gap_fracs[start:end]
            mean_H = sum(window_H) / window_size
            mean_gap = sum(window_gap) / window_size
            # use 1-based column indices
            out.write(f"{start+1}\t{end}\t{mean_H:.6f}\t{mean_gap:.6f}\n")

    print(f"[INFO] Window entropy saved to {win_tsv}", file=sys.stderr)
    print("[INFO] Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
