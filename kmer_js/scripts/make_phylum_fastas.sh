#!/bin/bash
set -e

# Phyla to extract (add or remove as needed)
PHYLA=(
  Arthropoda
  Chordata
  Mollusca
  Annelida
  Platyhelminthes
)

FASTA_IN="seqs_clean.fasta"

for P in "${PHYLA[@]}"; do
    echo "=== Phylum: $P ==="
    ./obtain_p_ids.sh "$P"
    IDS_FILE="${P}_ids.txt"
    OUT_FASTA="${P}.fasta"
    echo "Extracting $P sequences to $OUT_FASTA..."
    python extract_fasta_by_ids.py "$IDS_FILE" "$FASTA_IN" "$OUT_FASTA"
done

echo "Done. Per-phylum FASTA files generated."
