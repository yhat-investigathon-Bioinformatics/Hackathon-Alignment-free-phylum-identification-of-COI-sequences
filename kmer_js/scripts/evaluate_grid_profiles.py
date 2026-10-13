#!/usr/bin/env python3
# Usage:
#   python3 evaluate_grid_profiles.py > grid_eval.tsv
#
# Reads every profiles_*.json in ./profiles
# and computes:
#  - JS between each pair of phyla
#  - min_JS and mean_JS per (k,w,s) combination

import sys
import json
import glob
import math
from collections import defaultdict
from itertools import combinations

# Fixed list of phyla used in the grid
PHYLUM_LIST = [
    "Arthropoda",
    "Chordata",
    "Mollusca",
    "Annelida",
    "Platyhelminthes"
]

def js_divergence(p, q):
    """
    Jensen-Shannon divergence (in bits) between two distributions p and q.
    p and q: dicts {kmer: freq} with the same key support.
    """
    # Mixture M = (p + q) / 2
    m = {}
    for k in p.keys():
        m[k] = 0.5 * (p[k] + q[k])

    def kl(a, b):
        s = 0.0
        for k, v in a.items():
            if v <= 0.0:
                continue
            pb = b.get(k, 1e-12)  # smoothing
            s += v * math.log2(v / pb)
        return s

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

def parse_filename(fname):
    """
    Expects a name such as:
      profiles/profiles_Arthropoda_k4_w150_s10.json
    Returns: (phylum, k, w, s)
    """
    base = fname.split("/")[-1]       # profiles_Arthropoda_k4_w150_s10.json
    base = base.replace(".json", "")  # profiles_Arthropoda_k4_w150_s10
    parts = base.split("_")           # [profiles, Arthropoda, k4, w150, s10]

    phylum = parts[1]
    k = int(parts[2][1:])   # "k4" -> 4
    w = int(parts[3][1:])   # "w150" -> 150
    s = int(parts[4][1:])   # "s10" -> 10
    return phylum, k, w, s

def main():
    # Directory containing the profiles
    profiles_dir = "profiles"

    files = glob.glob(f"{profiles_dir}/profiles_*.json")
    if not files:
        print(f"No profiles_*.json files found in {profiles_dir}/", file=sys.stderr)
        sys.exit(1)

    # combos[(k,w,s)] = {phylum: profile_dict}
    combos = defaultdict(dict)

    for f in files:
        phylum, k, w, s = parse_filename(f)
        with open(f) as fh:
            prof = json.load(fh)
        combos[(k, w, s)][phylum] = prof

    # TSV header
    header = ["k", "window", "step", "min_JS", "mean_JS"]
    pair_names = []
    for a, b in combinations(PHYLUM_LIST, 2):
        pair_names.append(f"JS_{a}_{b}")
    header.extend(pair_names)

    print("\t".join(header))

    # Iterate over each combination
    for (k, w, s), phylum_profiles in sorted(combos.items()):
        # skip incomplete combinations
        if any(p not in phylum_profiles for p in PHYLUM_LIST):
            continue

        js_values = {}
        all_vals = []

        for a, b in combinations(PHYLUM_LIST, 2):
            pa = phylum_profiles[a]
            pb = phylum_profiles[b]
            js = js_divergence(pa, pb)
            js_values[(a, b)] = js
            all_vals.append(js)

        if not all_vals:
            continue

        min_js = min(all_vals)
        mean_js = sum(all_vals) / len(all_vals)

        row = [str(k), str(w), str(s), f"{min_js:.6f}", f"{mean_js:.6f}"]
        for a, b in combinations(PHYLUM_LIST, 2):
            row.append(f"{js_values[(a, b)]:.6f}")

        print("\t".join(row))

if __name__ == "__main__":
    main()
