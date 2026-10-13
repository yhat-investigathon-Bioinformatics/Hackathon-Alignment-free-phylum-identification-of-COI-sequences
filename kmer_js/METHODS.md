# Alignment-free phylum classification with k-mers and sliding windows

Methods and results of the k-mer + Jensen-Shannon approach.

---

## 1. Goal

Build a fully alignment-free phylum classifier for COI sequences based on:

- k-mer frequencies
- Sliding windows
- Jensen-Shannon (JS) divergence

The classifier should:

- Classify COI sequences of variable length
- Generalize to new sequences without prior alignment
- Be robust to fragmented or partial sequences

---

## 2. Dataset preparation

Inputs:

- `seqs_clean.fasta`: full FASTA, restricted to sequences containing only A/C/G/T
- `taxon.taxon`: taxonomic assignments

Sequences were split by phylum by matching the `p__<Phylum>` field, producing:

- `Arthropoda.fasta`
- `Chordata.fasta`
- `Mollusca.fasta`
- `Annelida.fasta`
- `Platyhelminthes.fasta`

---

## 3. Balanced subsampling

Arthropoda has more than 1,100,000 sequences, while the other phyla have 20,000–30,000.
To avoid bias and excessive compute cost:

`N_train = min(30,000, N_total)`

| Phylum | Sequences |
|--------|-----------|
| Arthropoda | 30,000 |
| Chordata | 30,000 |
| Mollusca | 30,000 |
| Annelida | 26,289 |
| Platyhelminthes | 24,911 |

---

## 4. Alignment-free approach

No alignments are used for classification. The classifier relies only on:

- k-mer distributions
- Relative frequencies
- Information-theoretic distances

This allows classifying short, long and incomplete sequences.

---

## 5. Sliding-window profiles

Each sequence is scanned with:

- k-mers of size **k**
- Windows of size **w**
- Step **s**

Each window yields a local k-mer profile. All local profiles are averaged into a global profile per
phylum, a vector of dimension 4^k (4,096 for k=6).

---

## 6. Parameter grid search

- k ∈ {4, 5, 6}
- w ∈ {150, 200, 250}
- s ∈ {10, 20}

18 combinations × 5 phyla = 90 profiles (≈ 8 h of compute).

---

## 7. Jensen-Shannon evaluation

For each combination, a 5×5 matrix of pairwise JS divergences between phyla was computed, and two
summary metrics were extracted:

- `min_JS`: minimum separation between any two phyla
- `mean_JS`: average separation

Results: [`results/grid_eval.tsv`](results/grid_eval.tsv).

![JS matrix k=6](../figures/kmer_js_matrix_k6.png)

---

## 8. Selecting the optimal setup

The goal is to maximize `min_JS`, so that even the closest pair of phyla stays well separated.

| k | min_JS |
|---|--------|
| 6 | ≈ 0.088 |
| 5 | ≈ 0.053 |

**k is the dominant parameter**: k=6 maximizes the separation between phyla. Among window and step
values, the selected setup is:

- k = 6
- w = 250
- s = 20

It offers near-maximal separation, lower compute cost and higher statistical stability.

---

## 9. PCA

PCA over the average profile of each phylum (k=6, w=250, s=20); each point is the compositional
signature of a phylum (`scripts/plot_profile_pca.py`).

- All phyla are separated
- Platyhelminthes is completely isolated
- Chordata is clearly separated from the protostome group
- Mollusca, Annelida and Arthropoda are closer to each other

---

## 10. JS distance matrix and clustering

Heatmap + Ward dendrogram of the JS distances between final profiles
(`scripts/plot_clustering_and_heatmap.py`).

- Smallest distances: Arthropoda–Annelida ≈ 0.088, Mollusca–Arthropoda ≈ 0.091
- Largest distance: Platyhelminthes–Chordata ≈ 0.341

Observed topology: Mollusca, Annelida and Arthropoda form a close block; Chordata separates from the
protostomes; Platyhelminthes is the most distant phylum.

---

## 11. Biological consistency

The pattern agrees with modern phylogeny:

- Lophotrochozoa: Mollusca + Annelida
- Ecdysozoa: Arthropoda
- Deuterostomia: Chordata
- Platyhelminthes: most divergent group

Mitochondrial genomes carry lineage-specific compositional biases, which explains the strong
separation obtained with k-mers.

---

## 12. Validation with a classifier

New profiles (k=6, w=250, s=20) were built from the training split of each phylum and evaluated with
`scripts/validate_classifier.py` on a held-out test set (200 sequences per phylum).

| True \ Predicted | Arthropoda | Chordata | Mollusca | Annelida | Platyhelminthes |
|---|---|---|---|---|---|
| **Arthropoda** | 187 | 2 | 3 | 2 | 4 |
| **Chordata** | 7 | 188 | 0 | 0 | 1 |
| **Mollusca** | 22 | 0 | 163 | 1 | 14 |
| **Annelida** | 32 | 3 | 11 | 152 | 1 |
| **Platyhelminthes** | 4 | 0 | 1 | 0 | 194 |

**Global accuracy: 89.1%** (5 phyla).

- Platyhelminthes (97%), Chordata (94%) and Arthropoda (93.5%) are classified almost perfectly.
- Mollusca (81.5%) is mostly confused with Arthropoda and Platyhelminthes.
- Annelida (76%) is the hardest class, mostly confused with Arthropoda, consistent with their small JS
  distance.

On the shared 8-phylum dataset used for the final comparison, this approach reaches **88.1%**
(`results/confusion_profiles_dl_ml_k6_w250_s20.json`).

![Confusion matrix, 8 phyla](../figures/kmer_confusion_8phyla.png)

---

## 13. Conclusion

The alignment-free system separates phyla using only k-mers, needs no alignment and works with
sequences of variable length. Final setup: **k = 6, w = 250, s = 20**.
