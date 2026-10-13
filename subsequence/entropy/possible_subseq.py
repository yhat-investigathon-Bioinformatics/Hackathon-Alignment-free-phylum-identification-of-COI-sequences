import pandas as pd

df = pd.read_csv("conservation_softtrim.tsv", sep="\t")
# columns: pos, coverage, conservation

window = 50

candidates = []
for start in range(1, df["pos"].max() - window + 2):
    end = start + window - 1
    sub = df[(df["pos"] >= start) & (df["pos"] <= end)]
    if (sub["coverage"] >= 0.85).mean() < 0.6:
        # require at least 60% of the window's columns
        # to have good coverage
        continue
    mean_cons = sub["conservation"].mean()
    candidates.append((start, end, mean_cons))

# sort by mean conservation (ascending: most variable first)
candidates_sorted = sorted(candidates, key=lambda x: x[2])

print("Most variable windows (low mean conservation):")
for c in candidates_sorted[:50]:
    print(f"{c[0]}-{c[1]}  mean_cons={c[2]:.3f}")

print("\nMost conserved windows (high mean conservation):")
for c in candidates_sorted[-5:]:
    print(f"{c[0]}-{c[1]}  mean_cons={c[2]:.3f}")
