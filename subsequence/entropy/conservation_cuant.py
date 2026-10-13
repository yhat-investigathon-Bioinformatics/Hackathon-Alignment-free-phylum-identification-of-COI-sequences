from Bio import SeqIO
import numpy as np

aln_file = "sample10k.softtrim.fasta"
out_tsv = "conservation_softtrim.tsv"

records = list(SeqIO.parse(aln_file, "fasta"))
n_seqs = len(records)
aln_len = len(records[0].seq)

print(f"Sequences: {n_seqs}, alignment length: {aln_len}")

arr = np.array([list(str(rec.seq)) for rec in records])

with open(out_tsv, "w") as out:
    out.write("pos\tcoverage\tconservation\n")
    for i in range(aln_len):
        col = arr[:, i]
        col_nogap = col[col != "-"]
        if len(col_nogap) == 0:
            coverage = 0.0
            conservation = 0.0
        else:
            coverage = len(col_nogap) / n_seqs
            unique, counts = np.unique(col_nogap, return_counts=True)
            conservation = counts.max() / len(col_nogap)
        out.write(f"{i+1}\t{coverage:.4f}\t{conservation:.4f}\n")

print(f"Saved to {out_tsv}")
