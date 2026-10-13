#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("results/grid_eval.tsv", sep="\t")

df["combo"] = (
    "k=" + df["k"].astype(str) +
    " w=" + df["window"].astype(str) +
    " s=" + df["step"].astype(str)
)

plt.figure(figsize=(12, 4))
sns.heatmap(
    df[["min_JS"]].T,
    cmap="viridis",
    cbar=True,
    yticklabels=["min_JS"],
    xticklabels=df["combo"],
)

plt.xticks(rotation=90)
plt.title("Minimum separation between phyla per combination")
plt.tight_layout()
plt.show()
