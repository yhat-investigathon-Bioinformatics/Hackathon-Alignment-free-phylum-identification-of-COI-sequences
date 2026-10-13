import json
import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# ===============================
# JENSEN-SHANNON
# ===============================

def js_divergence(p, q):
    all_keys = set(p.keys()) | set(q.keys())
    m = {k: 0.5 * (p.get(k, 0.0) + q.get(k, 0.0)) for k in all_keys}

    def kl(a, b):
        s = 0.0
        for k, v in a.items():
            if v <= 0:
                continue
            pb = b.get(k, 1e-12)
            s += v * math.log2(v / pb)
        return s

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

# ===============================
# k-MERS
# ===============================

def kmer_profile(seq, k):
    prof = {}
    total = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        prof[kmer] = prof.get(kmer, 0) + 1
        total += 1
    if total == 0:
        return {}
    for kmer in prof:
        prof[kmer] /= total
    return prof

def sliding_profile(seq, k, w, s):
    global_prof = {}
    n_windows = 0
    for i in range(0, len(seq) - w + 1, s):
        window = seq[i:i+w]
        local = kmer_profile(window, k)
        for kmer, val in local.items():
            global_prof[kmer] = global_prof.get(kmer, 0) + val
        n_windows += 1
    if n_windows == 0:
        return {}
    for kmer in global_prof:
        global_prof[kmer] /= n_windows
    return global_prof

# ===============================
# READ MULTI-SEQUENCE FASTA
# ===============================

def read_fasta_multi(fname):
    seqs = {}
    current_id = None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                seqs[current_id] = ""
            else:
                if current_id is not None:
                    seqs[current_id] += line.upper()
    return seqs

# ===============================
# DETECT PROFILES AND PARAMETERS (k, w, s)
# ===============================

def detect_phyla_and_params(profiles_dir):
    """
    Looks for files such as:
        ml_Annelida_k6_w250_s20.json
        ml_Arthropoda_k6_w250_s20.json
    and extracts:
        phylum, k, w, s

    Raises an error if k, w or s differ between files.
    """
    phyla = []
    params = set()
    files = []

    for fn in os.listdir(profiles_dir):
        if not fn.endswith(".json"):
            continue

        name = fn[:-5]  # strip .json
        parts = name.split("_")
        # expected form: [ml, Arthropoda, k6, w250, s20]
        if len(parts) != 5:
            continue

        prefix, phylum, k_part, w_part, s_part = parts
        if not (k_part.startswith("k") and w_part.startswith("w") and s_part.startswith("s")):
            continue

        try:
            k_val = int(k_part[1:])
            w_val = int(w_part[1:])
            s_val = int(s_part[1:])
        except ValueError:
            continue

        phyla.append(phylum)
        params.add((k_val, w_val, s_val))
        files.append((phylum, k_val, w_val, s_val, fn))

    if not files:
        raise RuntimeError(
            f"No valid profiles found in {profiles_dir} "
            "(expected names such as ml_Phylum_k6_w250_s20.json)"
        )

    if len(params) != 1:
        msg = "Error: multiple (k, w, s) combinations found among the profiles:\n"
        for (k_val, w_val, s_val) in sorted(params):
            msg += f"  k={k_val}, w={w_val}, s={s_val}\n"
        raise RuntimeError(msg)

    (k_val, w_val, s_val) = params.pop()
    phyla = sorted(set(phyla))

    return phyla, k_val, w_val, s_val, files

# ===============================
# MAIN
# ===============================

def main():
    BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    RESULTS_DIR = os.path.join(BASE_DIR, "results")
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if not (3 <= len(sys.argv) <= 4):
        print("Usage:")
        print("  python scripts/validate_classifier.py test_set.fasta test_labels.tsv [profiles_dir]")
        print("")
        print("  Examples:")
        print("    python scripts/validate_classifier.py data/splits/test_set.fasta data/splits/test_labels.tsv")
        print("    python scripts/validate_classifier.py data/splits/test_set.fasta data/splits/test_labels.tsv profiles/profiles_dl_ml")
        sys.exit(1)

    fasta_file = sys.argv[1]
    label_file = sys.argv[2]

    if len(sys.argv) == 4:
        profiles_dir_name = sys.argv[3]
    else:
        profiles_dir_name = os.path.join("profiles", "profiles_training")

    # Resolve profiles path
    if os.path.isabs(profiles_dir_name):
        PROFILES_DIR = profiles_dir_name
    else:
        PROFILES_DIR = os.path.join(BASE_DIR, profiles_dir_name)

    if not os.path.isdir(PROFILES_DIR):
        print(f"Error: profiles directory does not exist: {PROFILES_DIR}")
        sys.exit(1)

    print(f"Using profiles directory: {PROFILES_DIR}")

    # Detect phyla and k, w, s parameters from file names
    try:
        PHYLUMS, k, w, s, profile_files = detect_phyla_and_params(PROFILES_DIR)
    except RuntimeError as e:
        print(e)
        sys.exit(1)

    print("Detected phyla:", ", ".join(PHYLUMS))
    print(f"Detected parameters: k={k}, w={w}, s={s}")

    # Load trained profiles
    profiles = {}
    for ph in PHYLUMS:
        candidates = [fn for (ph2, k2, w2, s2, fn) in profile_files if ph2 == ph]
        if not candidates:
            print(f"Error: no profile file found for phylum {ph}")
            sys.exit(1)
        fname = os.path.join(PROFILES_DIR, candidates[0])
        with open(fname) as f:
            profiles[ph] = json.load(f)

    # Load test sequences
    seqs = read_fasta_multi(fasta_file)

    # Load true labels
    true_labels = {}
    with open(label_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            sid, lab = line.split()
            true_labels[sid] = lab

    # Confusion matrix
    confusion = {p: {q: 0 for q in PHYLUMS} for p in PHYLUMS}
    correct = 0
    total = 0

    # Ensure k-mer support (once)
    all_kmers = list(next(iter(profiles.values())).keys())

    for sid, seq in seqs.items():
        real = true_labels.get(sid)
        if real is None:
            continue

        if len(seq) < w:
            continue

        test_profile = sliding_profile(seq, k, w, s)
        if not test_profile:
            continue

        for kmer in all_kmers:
            test_profile.setdefault(kmer, 0.0)

        results = {}
        for ph in PHYLUMS:
            results[ph] = js_divergence(test_profile, profiles[ph])

        predicted = min(results, key=results.get)

        if real not in confusion:
            continue

        confusion[real][predicted] += 1
        total += 1
        if predicted == real:
            correct += 1

    acc = correct / total if total > 0 else 0.0

    # ===============================
    # TEXT REPORT
    # ===============================
    print("\nCONFUSION MATRIX:\n")
    print("\t" + "\t".join(PHYLUMS))
    for r in PHYLUMS:
        row = "\t".join(str(confusion[r][c]) for c in PHYLUMS)
        print(f"{r}\t{row}")
    print("\nOVERALL ACCURACY:", acc)

    # ===============================
    # SAVE JSON
    # ===============================
    tag = os.path.basename(os.path.normpath(PROFILES_DIR))
    json_path = os.path.join(RESULTS_DIR, f"confusion_{tag}_k{k}_w{w}_s{s}.json")

    # convert confusion dict to a matrix and save metadata
    mat = np.array([[confusion[r][c] for c in PHYLUMS] for r in PHYLUMS], dtype=int)

    out_data = {
        "phyla": PHYLUMS,
        "confusion_matrix": mat.tolist(),
        "confusion_dict": confusion,
        "k": k,
        "w": w,
        "s": s,
        "profiles_dir": PROFILES_DIR,
        "accuracy": acc,
        "normalized_by": 2000
    }

    with open(json_path, "w") as f:
        json.dump(out_data, f, indent=2)

    print(f"\nConfusion matrix JSON saved to: {json_path}")

    # ===============================
    # HEATMAP NORMALIZED BY 2000
    # ===============================
    mat_norm = mat.astype(float) / 2000.0

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(mat_norm)

    ax.set_xticks(range(len(PHYLUMS)))
    ax.set_yticks(range(len(PHYLUMS)))
    ax.set_xticklabels(PHYLUMS, rotation=45, ha="right")
    ax.set_yticklabels(PHYLUMS)

    for i in range(len(PHYLUMS)):
        for j in range(len(PHYLUMS)):
            ax.text(j, i, f"{mat_norm[i, j]:.3f}",
                    ha="center", va="center", fontsize=8)

    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Confusion matrix (normalized by 2000)")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()

    heatmap_path = os.path.join(RESULTS_DIR, f"confusion_heatmap_{tag}_k{k}_w{w}_s{s}.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.show()

    print(f"Heatmap saved to: {heatmap_path}")

if __name__ == "__main__":
    main()
