import pandas as pd
from Bio import SeqIO
import numpy as np
import random
import re

# ---------- CONFIG ----------
aln_file = "sample10k.softtrim.fasta"
tax_file = "taxon.taxon"
min_seqs = 20          # kept low so that all phyla are processed
max_seqs = 1000

# windows to evaluate (1-based)
windows = [
    (147, 196),
    (260, 330),
    (340, 390),
    (450, 530),
]

out_tsv = "per_phylum_window_stats.tsv"
# -----------------------------------

print("Reading alignment...")
records = list(SeqIO.parse(aln_file, "fasta"))
seqs = {rec.id: np.array(list(str(rec.seq))) for rec in records}
ids_aln = set(seqs.keys())
aln_len = len(records[0].seq)
print(f"MSA: {len(records)} sequences, {aln_len} columns")

print("Reading taxonomy...")
tax = pd.read_csv(tax_file, sep="\t", header=None, names=["seqid", "taxonomy"])

def extract_phylum(s):
    m = re.search(r"p__([A-Za-z0-9]+)", s)
    return m.group(1) if m else None

tax["phylum"] = tax["taxonomy"].apply(extract_phylum)

# keep only IDs present in the MSA
tax = tax[tax["seqid"].isin(ids_aln)]
print("Valid rows:", len(tax))

rows = []

# --- GROUP BY PHYLUM ---
for phylum, df_p in tax.groupby("phylum"):
    ids = df_p["seqid"].tolist()

    if len(ids) < min_seqs:
        continue

    # subsampling
    if len(ids) > max_seqs:
        ids = random.sample(ids, max_seqs)

    print(f"Processing phylum {phylum} ({len(ids)} seqs)...")

    mat = np.stack([seqs[i] for i in ids])  # shape (n, 680)

    for (start, end) in windows:
        cols = np.arange(start - 1, end)  # 0-based
        sub = mat[:, cols]

        # per-column coverage (non-gap)
        coverage = (sub != "-").mean(axis=0)

        # per-column conservation
        cons_cols = []
        for j in range(sub.shape[1]):
            col = sub[:, j]
            col_ng = col[col != "-"]
            if len(col_ng) == 0:
                continue
            vals, counts = np.unique(col_ng, return_counts=True)
            cons_cols.append(counts.max() / len(col_ng))

        if cons_cols:
            mean_cons = np.mean(cons_cols)
            frac_var = np.mean(np.array(cons_cols) < 0.7)
        else:
            mean_cons = np.nan
            frac_var = np.nan

        rows.append({
            "phylum": phylum,
            "subseq": f"{start}-{end}",
            "n_seqs": len(ids),
            "mean_conservation": mean_cons,
            "frac_columns_cons_below_0.7": frac_var,
            "mean_coverage": np.mean(coverage),
        })

# --- EXPORT RESULTS ---
df_out = pd.DataFrame(rows)
df_out.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved to {out_tsv} ({len(df_out)} rows)")
