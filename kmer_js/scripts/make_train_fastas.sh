#!/bin/bash

NMAX=30000

echo "[INFO] Generating training FASTA files with NMAX=${NMAX}"

# Arthropoda (very large)
seqkit sample -n ${NMAX} Arthropoda.fasta > Arthropoda_train.fasta

# Chordata (very large)
seqkit sample -n ${NMAX} Chordata.fasta > Chordata_train.fasta

# Mollusca (medium)
seqkit sample -n ${NMAX} Mollusca.fasta > Mollusca_train.fasta

# Annelida (already < 30k)
cp Annelida.fasta Annelida_train.fasta

# Platyhelminthes (already < 30k)
cp Platyhelminthes.fasta Platyhelminthes_train.fasta

echo "[DONE] Training FASTA files generated:"
ls *_train.fasta
