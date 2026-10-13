import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

phylums = [
    "Arthropoda",
    "Chordata",
    "Mollusca",
    "Annelida",
    "Platyhelminthes"
]

k, w, s = 6, 250, 20

X = []
labels = []

for ph in phylums:
    fname = f"profiles/profiles_{ph}_k{k}_w{w}_s{s}.json"
    with open(fname) as f:
        prof = json.load(f)
    X.append(list(prof.values()))
    labels.append(ph)

X = np.array(X)

pca = PCA(n_components=2)
X2 = pca.fit_transform(X)

plt.figure(figsize=(6,5))
for i, ph in enumerate(labels):
    plt.scatter(X2[i,0], X2[i,1], label=ph, s=80)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of phylum profiles (k=6, w=250, s=20)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
