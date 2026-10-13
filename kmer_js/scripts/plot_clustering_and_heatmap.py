import json
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# ===============================
# FIXED FINAL CONFIGURATION
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
PROFILES_DIR = "profiles"

# ===============================
# JENSEN-SHANNON DIVERGENCE
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
# LOAD PROFILES
# ===============================
profiles = {}

for ph in PHYLUMS:
    fname = f"{PROFILES_DIR}/profiles_{ph}_k{k}_w{w}_s{s}.json"
    with open(fname) as f:
        profiles[ph] = json.load(f)

# ===============================
# JS DISTANCE MATRIX
# ===============================
n = len(PHYLUMS)
dist_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        if i == j:
            dist_matrix[i, j] = 0.0
        else:
            dist_matrix[i, j] = js_divergence(
                profiles[PHYLUMS[i]],
                profiles[PHYLUMS[j]]
            )

# ===============================
# HIERARCHICAL CLUSTERING (WARD)
# ===============================
# Convert to condensed form for scipy
condensed_dist = squareform(dist_matrix)
Z = linkage(condensed_dist, method="ward")

# ===============================
# FIGURE 1: HEATMAP + DENDROGRAM
# ===============================
sns.set(style="white")

g = sns.clustermap(
    dist_matrix,
    row_linkage=Z,
    col_linkage=Z,
    xticklabels=PHYLUMS,
    yticklabels=PHYLUMS,
    cmap="viridis",
    linewidths=0.8,
    annot=True,
    fmt=".3f",
    figsize=(8, 8)
)

plt.suptitle(
    "Heatmap + hierarchical clustering (Ward)\n"
    f"JS divergence - k={k}, w={w}, s={s}",
    y=1.02
)

plt.show()

# ===============================
# FIGURE 2: WARD DENDROGRAM ONLY
# ===============================
plt.figure(figsize=(7, 4))

dendrogram(
    Z,
    labels=PHYLUMS,
    leaf_rotation=45,
    leaf_font_size=11,
    color_threshold=None
)

plt.title(f"Ward dendrogram (JS divergence)\nk={k}, w={w}, s={s}")
plt.ylabel("Distance")
plt.tight_layout()
plt.show()

print("Clustering and heatmap completed.")
