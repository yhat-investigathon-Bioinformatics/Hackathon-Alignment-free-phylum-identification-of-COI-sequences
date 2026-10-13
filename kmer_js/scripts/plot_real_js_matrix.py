import json
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
PROFILES_DIR = "profiles"

# ===============================
# JS FUNCTION
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
# LOAD OBSERVED PROFILES
# ===============================
profiles = {}

for ph in PHYLUMS:
    fname = f"{PROFILES_DIR}/profiles_{ph}_k{k}_w{w}_s{s}.json"
    with open(fname) as f:
        profiles[ph] = json.load(f)

# ===============================
# OBSERVED 5x5 JS MATRIX
# ===============================
n = len(PHYLUMS)
JS = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        if i == j:
            JS[i, j] = 0.0
        else:
            JS[i, j] = js_divergence(
                profiles[PHYLUMS[i]],
                profiles[PHYLUMS[j]]
            )

# ===============================
# PRINT NUMERIC MATRIX
# ===============================
print("\nOBSERVED JS DISTANCE MATRIX (k=6, w=250, s=20)\n")

print("\t" + "\t".join(PHYLUMS))
for i, ph in enumerate(PHYLUMS):
    values = "\t".join(f"{JS[i,j]:.6f}" for j in range(n))
    print(f"{ph}\t{values}")

# ===============================
# OBSERVED min_JS AND mean_JS
# ===============================
non_diag = JS[np.triu_indices(n, k=1)]

print("\nObserved min_JS:", non_diag.min())
print("Observed mean_JS:", non_diag.mean())

# ===============================
# HEATMAP
# ===============================
plt.figure(figsize=(7, 6))

sns.heatmap(
    JS,
    xticklabels=PHYLUMS,
    yticklabels=PHYLUMS,
    annot=True,
    fmt=".3f",
    cmap="viridis",
    linewidths=0.8
)

plt.title("JS distance matrix between phyla\nk=6, w=250, s=20")
plt.tight_layout()

# Save figure for the report
plt.savefig("js_matrix_real_k6.png", dpi=300)
plt.show()

print("\nFigure saved as: js_matrix_real_k6.png")
