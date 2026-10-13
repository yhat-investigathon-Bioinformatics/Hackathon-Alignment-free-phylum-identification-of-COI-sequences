#!/usr/bin/env python3
# ==========================================================
# build_profiles_sliding.py
#
# Usage:
#   python3 build_profiles_sliding.py input.fasta PHYLUM k window step out_prefix
#
# Example:
#   python3 build_profiles_sliding.py Arthropoda.fasta Arthropoda 5 200 20 profiles
#
# Output:
#   profiles_Arthropoda_k5_w200_s20.json
# ==========================================================
#
import sys
from collections import Counter
from itertools import product
import json

# ----------------------------
def all_kmers(k):
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=k)]

# ----------------------------
def sliding_kmer_freq(seq, k, window, step):
    """
    Return a list of dictionaries:
    [ {kmer: freq}, {kmer: freq}, ... ], one per valid window
    """
    seq = seq.upper().replace("N", "")
    L = len(seq)
    results = []

    if L < window or L < k:
        return results

    for start in range(0, L - window + 1, step):
        sub = seq[start:start + window]
        kmers = Counter()
        for i in range(len(sub) - k + 1):
            kmer = sub[i:i+k]
            if all(b in "ACGT" for b in kmer):
                kmers[kmer] += 1

        total = sum(kmers.values())
        if total > 0:
            freqs = {kmer: (cnt / total) for kmer, cnt in kmers.items()}
            results.append(freqs)

    return results

# ----------------------------
def fasta_iter(fname):
    with open(fname) as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, ''.join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, ''.join(seq)

# ----------------------------
def main():
    if len(sys.argv) != 7:
        print("Usage:")
        print("  python3 build_profiles_sliding.py input.fasta PHYLUM k window step out_prefix")
        sys.exit(1)

    fasta_file = sys.argv[1]
    phylum     = sys.argv[2]
    k          = int(sys.argv[3])
    window     = int(sys.argv[4])
    step       = int(sys.argv[5])
    out_prefix = sys.argv[6]

    allkm = all_kmers(k)

    global_sum = Counter()
    n_windows = 0
    n_seqs = 0

    print(f"[INFO] Phylum: {phylum}")
    print(f"[INFO] k={k}, window={window}, step={step}")
    print(f"[INFO] Processing {fasta_file} ...")

    for header, seq in fasta_iter(fasta_file):
        n_seqs += 1
        windows = sliding_kmer_freq(seq, k, window, step)
        for w in windows:
            for km, v in w.items():
                global_sum[km] += v
            n_windows += 1

        if n_seqs % 10000 == 0:
            print(f"[INFO] {n_seqs} sequences processed...")

    if n_windows == 0:
        print("[ERROR] No valid windows were generated.")
        sys.exit(1)

    profile = {km: global_sum.get(km, 0) / n_windows for km in allkm}

    out_file = f"{out_prefix}_{phylum}_k{k}_w{window}_s{step}.json"

    with open(out_file, "w") as f:
        json.dump(profile, f, indent=2)

    print(f"[DONE] Profile saved to: {out_file}")
    print(f"[STATS] Sequences: {n_seqs}")
    print(f"[STATS] Total windows: {n_windows}")

if __name__ == "__main__":
    main()
