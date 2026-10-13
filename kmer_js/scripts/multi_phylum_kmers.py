#!/usr/bin/env python3
# Usage:
#   python3 multi_phylum_kmers.py k out_prefix \
#       Arthropoda.fasta:Arthropoda \
#       Chordata.fasta:Chordata \
#       Mollusca.fasta:Mollusca ...

import sys
import math
from collections import Counter
from itertools import product
import json

def all_kmers(k):
    bases = ['A','C','G','T']
    return [''.join(p) for p in product(bases, repeat=k)]

def seq_kmer_freqs(seq, k):
    seq = seq.upper().replace('N', '')
    L = len(seq)
    if L < k:
        return {}
    cnt = Counter(seq[i:i+k] for i in range(L-k+1))
    # remove k-mers with non-ACGT bases
    for km in list(cnt.keys()):
        if any(b not in 'ACGT' for b in km):
            del cnt[km]
    total = sum(cnt.values())
    if total == 0:
        return {}
    return {km: c / total for km, c in cnt.items()}

def aggregate_profile(fasta, k, max_seqs=200000):
    """
    Read a FASTA file and return:
      - mean k-mer profile (dict)
      - number of sequences used

    max_seqs: if not None, limits the number of sequences (useful for very large phyla).
    """
    totals = Counter()
    nseq = 0
    with open(fasta) as f:
        hdr = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # close previous sequence
                if hdr is not None:
                    s = "".join(seq_lines)
                    freqs = seq_kmer_freqs(s, k)
                    for km, v in freqs.items():
                        totals[km] += v
                    nseq += 1
                    if max_seqs is not None and nseq >= max_seqs:
                        break
                hdr = line
                seq_lines = []
            else:
                seq_lines.append(line)
        # last sequence, if max_seqs was not reached
        if hdr is not None and (max_seqs is None or nseq < max_seqs):
            s = "".join(seq_lines)
            freqs = seq_kmer_freqs(s, k)
            for km, v in freqs.items():
                totals[km] += v
            nseq += 1

    if nseq == 0:
        raise ValueError(f"No valid sequences in {fasta}")

    kmers = all_kmers(k)
    prof = {km: totals.get(km, 0.0) / nseq for km in kmers}
    return prof, nseq

def kl_div(a, b):
    s = 0.0
    for k, v in a.items():
        if v <= 0.0:
            continue
        pb = b.get(k, 1e-12)
        s += v * math.log2(v / pb)
    return s

def js_divergence(p, q):
    kmers = p.keys()
    m = {k: 0.5 * (p.get(k, 0.0) + q.get(k, 0.0)) for k in kmers}
    return 0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 multi_phylum_kmers.py k out_prefix fasta1:Name1 [fasta2:Name2 ...]")
        sys.exit(1)

    k = int(sys.argv[1])
    out_prefix = sys.argv[2]
    group_specs = sys.argv[3:]  # list of "path:Name"

    # Set this to limit the number of sequences per phylum,
    # e.g. 200000 to reduce runtime/RAM on Arthropoda
    MAX_SEQS_PER_GROUP = None  # or an int, e.g. 200000

    profiles = {}
    nseqs = {}

    # 1) Build profiles per group
    for spec in group_specs:
        path, name = spec.split(":", 1)
        print(f"[INFO] Processing {name} from {path} ...", file=sys.stderr)
        prof, n = aggregate_profile(path, k, max_seqs=MAX_SEQS_PER_GROUP)
        profiles[name] = prof
        nseqs[name] = n

    # 2) Summary
    print("Group summary:", file=sys.stderr)
    for name in profiles:
        print(f"  {name}: nseq={nseqs[name]}", file=sys.stderr)

    # 3) JS matrix between groups
    names = list(profiles.keys())
    js_matrix = {a: {} for a in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                js = 0.0
            else:
                js = js_divergence(profiles[a], profiles[b])
            js_matrix[a][b] = js

    # Save the full JS matrix as a long-format TSV
    with open(out_prefix + "_js_matrix.tsv", "w") as out:
        out.write("phylum1\tphylum2\tJS_bits\n")
        for a in names:
            for b in names:
                out.write(f"{a}\t{b}\t{js_matrix[a][b]:.6f}\n")

    # 4) Save full profiles and counts
    with open(out_prefix + "_profiles.json", "w") as f:
        json.dump({"k": k, "profiles": profiles, "nseqs": nseqs}, f)
