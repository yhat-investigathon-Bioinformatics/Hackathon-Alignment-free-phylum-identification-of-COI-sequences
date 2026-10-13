#!/bin/bash

# Usage: ./obtain_p_ids.sh <Phylum>
# Example: ./obtain_p_ids.sh Arthropoda

if [ -z "$1" ]; then
    echo "Error: a phylum name is required."
    echo "Usage: $0 <Phylum>"
    exit 1
fi

PHYLUM="$1"

# taxonomy file
TAXON_FILE="taxon.taxon"

# output file
OUT_FILE="${PHYLUM}_ids.txt"

echo "Searching for p__${PHYLUM} in ${TAXON_FILE}..."

grep "p__${PHYLUM}" "$TAXON_FILE" | cut -f1 > "$OUT_FILE"

N=$(wc -l < "$OUT_FILE")

echo "Done. Found $N sequences."
echo "Saved to: $OUT_FILE"
grep "p__x" taxon.taxon | cut -f1 > x_ids.txt
