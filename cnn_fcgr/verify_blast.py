import pandas as pd
import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# ==========================================
# 1. CONFIGURATION
# ==========================================
PATH_PREDICTIONS = "submission_final_smart.csv"
PATH_QUERY_FASTA = "query.fasta"
SAMPLES_PER_PHYLUM = 3  # Number of sequences to check per phylum

def main():
    print("=== BLAST VERIFICATION (STRATIFIED BY PHYLUM) ===")

    # 1. Load data
    print("Loading predictions...")
    df = pd.read_csv(PATH_PREDICTIONS)

    # Drop "Unclassified" / "Rejected" entries
    df_clean = df[~df['Final_Prediction'].str.contains("Unclassified|Rejected|Invalid")]

    # 2. Load FASTA sequences
    print("Indexing FASTA sequences...")
    seq_dict = {rec.id.strip(): str(rec.seq) for rec in SeqIO.parse(PATH_QUERY_FASTA, "fasta")}

    # 3. STRATIFIED SELECTION
    # Take 'n' samples from EACH available class
    target_ids = []

    # Group by predicted class
    groups = df_clean.groupby("Final_Prediction")

    print(f"\nSampling strategy: {SAMPLES_PER_PHYLUM} sequences per phylum.")

    for phylum, group in groups:
        # Sort by confidence (test the most confident predictions)
        # For random sampling, replace with .sample()
        top_samples = group.sort_values("Confidence", ascending=False).head(SAMPLES_PER_PHYLUM)

        group_ids = top_samples["Sequence_ID"].tolist()
        target_ids.extend(group_ids)

        print(f"   -> {phylum}: {len(group_ids)} samples selected.")

    print(f"\nTOTAL TO VERIFY: {len(target_ids)} sequences.")
    print("   (Approx. 30-60 s per sequence)\n")

    # 4. BLAST loop
    matches = 0
    total_checked = 0

    for i, seq_id in enumerate(target_ids):
        row = df[df["Sequence_ID"] == seq_id].iloc[0]
        model_prediction = row["Final_Prediction"]
        model_confidence = row["Confidence"]
        sequence = seq_dict.get(seq_id)

        if not sequence:
            print(f"Error: ID {seq_id} not found in FASTA.")
            continue

        print(f"[{i+1}/{len(target_ids)}] ID: {seq_id}")
        print(f"   Model prediction: \033[94m{model_prediction}\033[0m (Conf: {model_confidence:.2f})") # Blue
        print("   Querying NCBI... ", end="", flush=True)

        try:
            # Submit to BLAST
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1)
            blast_record = NCBIXML.read(result_handle)

            if blast_record.alignments:
                top_hit = blast_record.alignments[0]
                hit_def = top_hit.title

                # Extract organism name for display
                organism = hit_def.split("|")[-1].strip()[:100]

                print(f"\nBLAST TOP HIT: \033[92m{organism}\033[0m") # Green

                # Simple case-insensitive match check
                if model_prediction.lower() in hit_def.lower():
                    print("MATCH")
                    matches += 1
                else:
                    print("   MISMATCH (check manually)")

            else:
                print("\n   No BLAST hits.")

        except Exception as e:
            print(f"\n   Connection error: {e}")

        print("-" * 60)
        total_checked += 1
        time.sleep(5)

if __name__ == "__main__":
    main()
