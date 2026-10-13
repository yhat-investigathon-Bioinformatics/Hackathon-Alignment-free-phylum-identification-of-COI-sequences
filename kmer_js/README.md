# k-mer + Jensen-Shannon classifier (alignment-free)

Each phylum is represented by the average k-mer frequency profile of its training sequences,
computed over sliding windows. A new sequence is assigned to the phylum whose profile has the lowest
Jensen-Shannon divergence to it.

- Methods and results: [`METHODS.md`](METHODS.md)
- Description of every script: [`SCRIPTS.md`](SCRIPTS.md)

## Layout

```
kmer_js/
├── scripts/                     # all scripts (see SCRIPTS.md)
├── profiles/
│   ├── profiles_training/       # final profiles, 5 phyla (k=6, w=250, s=20)
│   └── profiles_dl_ml/          # final profiles, shared 8-phylum dataset
└── results/
    ├── grid_eval.tsv            # grid search: k, window, step, min_JS, mean_JS
    └── confusion_profiles_dl_ml_k6_w250_s20.json
```

Input data (not versioned) is expected under `kmer_js/data/`:
`data/raw/` (MIDORI2 FASTA + taxonomy), `data/per_phylum/` and `data/splits/` (train/test FASTAs,
`test_set.fasta`, `test_labels.tsv`).

## Pipeline

```
MIDORI2 FASTA + taxonomy
        │
        ▼
Cleaning + filtering (length, non-ACGT)      obtain_p_ids.sh, make_phylum_fastas.sh
        │
        ▼
data/per_phylum/*.fasta
        │
        ▼
80/20 train/test split                       make_train_fastas.sh
        │
        ▼
Profile building (k=6, w=250, s=20)          build_profiles_sliding.py
        │
        ▼
profiles/profiles_training/*.json
        │
        ▼
Classification                               classify_sequence.py
        │
        ▼
Evaluation: confusion matrix + accuracy      validate_classifier.py
```

## Quick start

```bash
cd kmer_js
python scripts/classify_sequence.py query.fasta
```
