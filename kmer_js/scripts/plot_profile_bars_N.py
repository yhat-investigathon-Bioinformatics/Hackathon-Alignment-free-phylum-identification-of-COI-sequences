import json
import matplotlib.pyplot as plt

# =========================
# CONFIGURATION
# =========================

fname = "profiles/profiles_Arthropoda_k6_w250_s20.json"
N = 50   # number of k-mers to show (30, 50, 100, etc.)

# =========================
# LOAD PROFILE
# =========================

with open(fname) as f:
    prof = json.load(f)

# =========================
# TAKE TOP N BY FREQUENCY
# =========================

top_items = sorted(prof.items(), key=lambda x: -x[1])[:N]

# =========================
# SORT THOSE N ALPHABETICALLY (X AXIS)
# =========================

top_items_alpha = sorted(top_items, key=lambda x: x[0])

kmers = [x[0] for x in top_items_alpha]
vals  = [x[1] for x in top_items_alpha]

# =========================
# PLOT
# =========================

plt.figure(figsize=(max(12, N/4), 5))
plt.bar(kmers, vals)
plt.xticks(rotation=90)
plt.ylabel("Frequency")
plt.xlabel("k-mer (alphabetical order)")
plt.title(f"Top {N} k-mers - Arthropoda (k=6, sorted alphabetically)")
plt.tight_layout()
plt.show()
