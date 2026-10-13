#!/usr/bin/env python3
# ==========================================================
# grid_search_sliding.py (TRAINING VERSION)
# Uses the balanced *_train.fasta files
# ==========================================================

import subprocess

# ----------------------------
# CONFIGURATION
# ----------------------------

phylums = {
    "Arthropoda": "Arthropoda_train.fasta",
    "Chordata": "Chordata_train.fasta",
    "Mollusca": "Mollusca_train.fasta",
    "Annelida": "Annelida_train.fasta",
    "Platyhelminthes": "Platyhelminthes_train.fasta"
}

k_values      = [4, 5, 6]
window_values = [150, 200, 250]
step_values   = [10, 20]

out_prefix = "profiles"

# ----------------------------
print("[INFO] Starting grid search with training FASTA files")
print("[INFO] Phylums:", ", ".join(phylums.keys()))
print("[INFO] k:", k_values)
print("[INFO] window:", window_values)
print("[INFO] step:", step_values)
print("--------------------------------------------------")

# ----------------------------
for k in k_values:
    for window in window_values:
        for step in step_values:
            print(f"\n[GRID] k={k}, window={window}, step={step}")

            for phylum, fasta in phylums.items():

                cmd = [
                    "python3",
                    "build_profiles_sliding.py",
                    fasta,
                    phylum,
                    str(k),
                    str(window),
                    str(step),
                    out_prefix
                ]

                print("  ->", phylum)
                subprocess.run(cmd, check=True)

print("\n[DONE] Grid search with training FASTA files completed.")
