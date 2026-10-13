# Data

The data is not versioned because of its size.

1. Download the `MIDORI2_UNIQ_NUC_GB268_CO1` set in FASTA format (`.fasta` + `.taxon`) from MIDORI2:
   <https://www.reference-midori.info/download.php>
2. Copy both files into `cnn_fcgr/` and run `python preprocessing.py`. It generates
   `dataset_10k_{train,test}.{fasta,taxon}` (8 phyla × 10,000 sequences, 80/20 split, seed 42).
