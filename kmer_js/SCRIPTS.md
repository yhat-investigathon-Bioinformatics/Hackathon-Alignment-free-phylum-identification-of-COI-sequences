# Scripts overview

## Core model

| Script | Description |
|---|---|
| `build_profiles_sliding.py` | Builds a k-mer profile per phylum using sliding windows. Input: training FASTA. Output: JSON with normalized k-mer counts. |
| `classify_sequence.py` | Classifies one or more FASTA sequences. Returns the JS divergence against each phylum and the predicted phylum. |
| `validate_classifier.py` | Evaluates the classifier on the held-out test set. Output: confusion matrix + global accuracy. |

## Grid search

| Script | Description |
|---|---|
| `grid_search_sliding.py` | Runs all combinations of k (k-mer size), w (window) and s (step), building profiles and the JS matrix for each. |
| `evaluate_grid_profiles.py` | Computes `min_JS` and `mean_JS` for each parameter combination (`results/grid_eval.tsv`). |
| `heatmap_grid.py` | Heatmap of grid search performance. |

## Plots

| Script | Description |
|---|---|
| `plot_profile_pca.py` | PCA of the per-phylum profiles. |
| `plot_profile_bars.py` | Bar plot of the most frequent k-mers. |
| `plot_profile_bars_N.py` | Same, for a configurable top N. |
| `plot_real_js_matrix.py` | Computes and plots the JS matrix between phyla. |
| `plot_clustering_and_heatmap.py` | Hierarchical clustering (Ward linkage) + heatmap. |
| `lengths_graph.py` | Histogram of sequence lengths. |

## Preprocessing

| Script | Description |
|---|---|
| `obtain_p_ids.sh` | Extracts sequence IDs per phylum from `taxon.taxon`. |
| `make_phylum_fastas.sh` | Builds one FASTA per phylum from the clean dataset. |
| `make_train_fastas.sh` | 80/20 split into `*_train.fasta` and `*_test.fasta`. |
| `split_10k_by_phylum.py` | Splits the shared 10k-per-phylum dataset into per-phylum FASTAs. |
| `extract_fasta_by_ids.py` | Extracts sequences from a list of IDs. |
| `existing_phylums.py` | Lists the phyla present in a FASTA file. |
| `multi_phylum_kmers.py` | Builds profiles for several phyla in one run. |

## Entropy

| Script | Description |
|---|---|
| `shannon_entropy_msa.py` | Per-column Shannon entropy of an MSA. Output: TSV with entropy and % gaps. |
