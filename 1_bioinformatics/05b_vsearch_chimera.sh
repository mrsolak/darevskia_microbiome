#!/usr/bin/env bash

PROJECT_DIR="/mnt/d/Rwd_claude/darevskia_16S_first"

vsearch \
    --uchime_ref  "$PROJECT_DIR/phyloseqs/asvs_sanger.fasta" \
    --db          "$PROJECT_DIR/raw_files/silva_nr99_v138.1_train_set.fa.gz" \
    --nonchimeras "$PROJECT_DIR/phyloseqs/non_chimeras_sanger.fasta" \
    --chimeras    "$PROJECT_DIR/phyloseqs/chimeras_sanger.fasta" \
    --uchimeout   "$PROJECT_DIR/results/tables/uchime_sanger.txt" \
    --threads 8 \
    --strand plus

echo "Done. Non-chimeras -> phyloseqs/non_chimeras_sanger.fasta"
