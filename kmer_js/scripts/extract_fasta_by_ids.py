import sys
if len(sys.argv) != 4:
    print("Usage: python extract_fasta_by_ids.py ids.txt input.fasta output.fasta")
    sys.exit(1)

ids_file = sys.argv[1]
fasta_in = sys.argv[2]
fasta_out = sys.argv[3]

# Read IDs into a set
with open(ids_file) as f:
    target_ids = {line.strip() for line in f if line.strip()}

print(f"Target IDs: {len(target_ids)}")

with open(fasta_in) as fin, open(fasta_out, "w") as fout:
    keep = False
    current_id = None

    for line in fin:
        if line.startswith(">"):
            # ID = first word of the header, without '>'
            current_id = line[1:].strip().split()[0]
            keep = current_id in target_ids
        if keep:
            fout.write(line)
