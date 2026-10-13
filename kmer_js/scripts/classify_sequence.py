import sys
import json
import math

# ===============================
# FINAL CONFIGURATION
# ===============================

PHYLUMS = [
    "Arthropoda",
    "Chordata",
    "Mollusca",
    "Annelida",
    "Platyhelminthes"
]

k = 6
w = 250
s = 20
PROFILES_DIR = "profiles/profiles_training"

# ===============================
# JENSEN-SHANNON
# ===============================

def js_divergence(p, q):
    m = {k: 0.5 * (p[k] + q[k]) for k in p}

    def kl(a, b):
        s = 0.0
        for k, v in a.items():
            if v <= 0:
                continue
            pb = b.get(k, 1e-12)
            s += v * math.log2(v / pb)
        return s

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

# ===============================
# k-MER COUNTING
# ===============================

def kmer_profile(seq, k):
    prof = {}
    total = 0

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        prof[kmer] = prof.get(kmer, 0) + 1
        total += 1

    for kmer in prof:
        prof[kmer] /= total

    return prof

# ===============================
# SLIDING WINDOW
# ===============================

def sliding_profile(seq, k, w, s):
    global_prof = {}
    n_windows = 0

    for i in range(0, len(seq) - w + 1, s):
        window = seq[i:i+w]
        local = kmer_profile(window, k)

        for kmer, val in local.items():
            global_prof[kmer] = global_prof.get(kmer, 0) + val

        n_windows += 1

    for kmer in global_prof:
        global_prof[kmer] /= n_windows

    return global_prof

# ===============================
# READ MULTI-SEQUENCE FASTA
# ===============================

def read_fasta_multi(fname):
    seqs = {}
    current_id = None

    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                seqs[current_id] = ""
            else:
                seqs[current_id] += line.upper()

    return seqs

# ===============================
# MAIN
# ===============================

if len(sys.argv) != 2:
    print("Usage:")
    print("  python classify_sequence.py sequences.fasta")
    sys.exit(1)

fasta_file = sys.argv[1]

# Load trained profiles
profiles = {}

for ph in PHYLUMS:
    fname = f"{PROFILES_DIR}/profiles_{ph}_k{k}_w{w}_s{s}.json"
    with open(fname) as f:
        profiles[ph] = json.load(f)

# Read all sequences
seqs = read_fasta_multi(fasta_file)

print("\n========================================")
print("MULTI-SEQUENCE CLASSIFICATION")
print("========================================")

for seq_id, seq in seqs.items():

    print(f"\nSequence: {seq_id}")
    print(f"Length: {len(seq)} bp")

    if len(seq) < w:
        print("[WARN] Too short to classify with w=250")
        continue

    test_profile = sliding_profile(seq, k, w, s)

    # Ensure full k-mer support
    all_kmers = list(next(iter(profiles.values())).keys())
    for kmer in all_kmers:
        test_profile.setdefault(kmer, 0.0)

    results = {}

    for ph in PHYLUMS:
        js = js_divergence(test_profile, profiles[ph])
        results[ph] = js

    for ph in PHYLUMS:
        print(f"{ph:20s} JS = {results[ph]:.6f}")

    predicted = min(results, key=results.get)
    print("PREDICTED PHYLUM:", predicted)

print("\n========================================")
print("CLASSIFICATION FINISHED")
print("========================================")
