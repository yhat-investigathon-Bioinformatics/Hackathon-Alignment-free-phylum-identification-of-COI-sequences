import numpy as np
import pandas as pd
import re
import os
import tensorflow as tf
from Bio import SeqIO
from tqdm import tqdm
from sklearn.preprocessing import LabelEncoder

# ==========================================
# 1. CONFIGURATION
# ==========================================
# Input/output files
PATH_QUERY_FASTA = "query.fasta"   # Query file
PATH_MODEL = "classifier.keras"
PATH_CLASSES = "classes.npy"
OUTPUT_CSV = "submission_final_smart.csv"   # Final .csv

WINDOW_LEN = 150  # Window size
STRIDE = 50       # Scan step (overlap)
K_MER = 6         # 2^6 = 64x64 pixels (must match the model)

# DECISION THRESHOLD
CONFIDENCE_THRESHOLD = 0.75

# ==========================================
# 2. FCGR COMPUTATION
# ==========================================
def sequence_to_fcgr(sequence, k):
    """Converts a DNA sequence into an FCGR image tensor."""
    sequence = sequence.upper()
    s_clean = re.sub(r'[^ACGT]', '', sequence)

    # Too short (< k): return None
    if len(s_clean) < k: return None

    # Fast base-to-digit mapping
    table = str.maketrans("ACGT", "0123")
    try:
        s_nums = [int(c) for c in s_clean.translate(table)]
    except:
        return None

    dim = 2**k
    matrix = np.zeros((dim, dim), dtype=np.float32)
    powers = [2**(k - 1 - j) for j in range(k)]

    for i in range(len(s_nums) - k + 1):
        kmer_nums = s_nums[i:i+k]
        row, col = 0, 0
        for j, base in enumerate(kmer_nums):
            p = powers[j]
            if base == 2: row += p; col += p
            elif base == 3: row += p
            elif base == 1: col += p
        matrix[row, col] += 1

    matrix = np.log1p(matrix)
    if np.max(matrix) > 0: matrix = matrix / np.max(matrix)
    return matrix

def load_classes():
    """Restores the exact phylum order from classes.npy."""
    if not os.path.exists(PATH_CLASSES):
        if os.path.exists("y_train.npy"):
            print("classes.npy not found, loading from y_train.npy...")
            y_train_labels = np.load("y_train.npy")
            le = LabelEncoder()
            le.fit(y_train_labels)
            return le
        else:
            raise FileNotFoundError("classes.npy is missing: cannot map model outputs to class names.")

    # Load saved classes
    print("Loading classes from classes.npy...")
    class_names = np.load(PATH_CLASSES)
    le = LabelEncoder()
    le.fit(class_names) # Fit the encoder to the saved classes
    print(f"   Classes loaded: {len(le.classes_)}")
    return le

# ==========================================
# 3. SCANNING LOGIC
# ==========================================
def scan_and_predict_smart(full_seq, model):
    """
    Strategy: sliding-window scan to find the best-matching window.
    """
    length = len(full_seq)

    # --- STRATEGY A: FALLBACK (short sequence) ---
    if length < WINDOW_LEN:
        img = sequence_to_fcgr(full_seq, K_MER)
        if img is None:
            return None, 0.0, "Invalid_Length"

        # Single inference
        x = np.expand_dims(img, axis=0)      # Batch of 1
        x = np.expand_dims(x, axis=-1)       # Channel
        probs = model.predict(x, verbose=0)[0]
        return np.argmax(probs), np.max(probs), "Short_WholeSeq"

    # --- STRATEGY B: SCANNING (long sequence) ---
    window_imgs = []
    window_starts = []

    # Generate windows
    for i in range(0, length - WINDOW_LEN + 1, STRIDE):
        roi = full_seq[i : i + WINDOW_LEN]
        img = sequence_to_fcgr(roi, K_MER)
        if img is not None:
            window_imgs.append(img)
            window_starts.append(i)

    if not window_imgs:
        return None, 0.0, "Invalid_Content"

    # Batch inference
    X_batch = np.array(window_imgs, dtype=np.float32)
    X_batch = np.expand_dims(X_batch, axis=-1)

    all_probs = model.predict(X_batch, verbose=0)

    # Find the window with the highest confidence
    max_conf_per_window = np.max(all_probs, axis=1) # Max probability per window
    best_window_idx_arr = np.argmax(max_conf_per_window) # Array index

    best_confidence = max_conf_per_window[best_window_idx_arr]
    best_class_idx = np.argmax(all_probs[best_window_idx_arr])

    # Actual position in the sequence (for the report)
    real_pos_start = window_starts[best_window_idx_arr]

    return best_class_idx, best_confidence, f"Scan_Pos_{real_pos_start}"

# ==========================================
# 4. MAIN
# ==========================================
def main():
    print(f"=== STARTING PREDICTOR ===")
    print(f"Reading: {PATH_QUERY_FASTA}")
    print(f"Model: {PATH_MODEL}")
    print(f"Confidence threshold: {CONFIDENCE_THRESHOLD}")

    # 1. Load resources
    try:
        le = load_classes()
        model = tf.keras.models.load_model(PATH_MODEL)
        print("Model loaded successfully.")
    except Exception as e:
        print(f"Error loading resources: {e}")
        return

    # 2. Report lists
    ids, best_guesses, confidences, final_decisions, methods = [], [], [], [], []

    # 3. Process FASTA file
    # Count first for the progress bar
    records = list(SeqIO.parse(PATH_QUERY_FASTA, "fasta"))
    total_seqs = len(records)

    print(f"Processing {total_seqs} sequences...")

    for record in tqdm(records):
        seq_id = record.id.strip()
        full_seq = str(record.seq).upper()

        idx_pred, conf, method = scan_and_predict_smart(full_seq, model)

        if idx_pred is not None:
            # Map index to name (0 -> Arthropoda)
            phylum_name = le.inverse_transform([idx_pred])[0]

            ids.append(seq_id)
            best_guesses.append(phylum_name)
            confidences.append(conf)
            methods.append(method)

            if conf >= CONFIDENCE_THRESHOLD:
                final_decisions.append(phylum_name)
            else:
                final_decisions.append("Unclassified (Low Conf)")

        else:
            # Error case
            ids.append(seq_id)
            best_guesses.append("Invalid")
            confidences.append(0.0)
            final_decisions.append("Rejected")
            methods.append(method)

    # 4. Save results
    df = pd.DataFrame({
        "Sequence_ID": ids,
        "Final_Prediction": final_decisions,
        "Best_Guess": best_guesses,
        "Confidence": confidences,
        "Method_Info": methods
    })

    # Sort: most confident first
    df = df.sort_values(by="Confidence", ascending=False)

    df.to_csv(OUTPUT_CSV, index=False)

    print("\n" + "="*40)
    print(f"Report generated: {OUTPUT_CSV}")
    print("="*40)
    print("\nPrediction summary:")
    print(df['Final_Prediction'].value_counts().head(10))
    print(f"\nTotal processed: {total_seqs}")

if __name__ == "__main__":
    main()
