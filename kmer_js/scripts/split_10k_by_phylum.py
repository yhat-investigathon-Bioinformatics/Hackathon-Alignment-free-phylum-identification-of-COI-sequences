import sys
import os

if len(sys.argv) != 3:
    print("Usage:")
    print("  python scripts/split_10k_by_phylum.py dataset_10k_train.fasta dataset_10k_train.taxon")
    sys.exit(1)

fasta_path = sys.argv[1]
taxon_path = sys.argv[2]

# 1) Read the .taxon file -> dict: id -> phylum
id2phylum = {}

with open(taxon_path) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        seq_id = parts[0]       # sequence ID (same as in the FASTA)
        phylum = parts[1]       # phylum name (e.g. Arthropoda)
        id2phylum[seq_id] = phylum

# 2) Open output files on demand
handles = {}

def get_handle(phylum):
    if phylum not in handles:
        out_name = f"{phylum}_10k_train.fasta"
        handles[phylum] = open(out_name, "w")
    return handles[phylum]

# 3) Iterate over the mixed FASTA and write each sequence to its phylum file
with open(fasta_path) as f:
    current_handle = None

    for line in f:
        if line.startswith(">"):
            header = line[1:].strip()
            seq_id = header.split()[0]   # use the first word as the ID
            phylum = id2phylum.get(seq_id)

            if phylum is None:
                # ID not found in the .taxon file -> skip it
                current_handle = None
                continue

            current_handle = get_handle(phylum)
            current_handle.write(line)

        else:
            if current_handle is not None:
                current_handle.write(line)

# 4) Close all files
for h in handles.values():
    h.close()

print("Done. FASTA files split by phylum:")
for ph in sorted(handles.keys()):
    print(f"  {ph}_10k_train.fasta")
