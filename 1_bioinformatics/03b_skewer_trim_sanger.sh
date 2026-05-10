#!/usr/bin/env bash
# Trim 16S primers from demultiplexed reads using bare trus/trus_ sequences.
# -m head removes everything up to and including the primer from the 5' end,
# which also removes any inline barcode that sits before the primer.
# Input:  raw_files/1_demultiplexed_sanger/*-assigned-trustrus_-pair[12].fastq
# Output: raw_files/2_trimmed_sanger/

PROJECT_DIR="/mnt/d/Rwd_claude/darevskia_16S_first"
IN_DIR="$PROJECT_DIR/raw_files/1_demultiplexed_sanger"
OUT_DIR="$PROJECT_DIR/raw_files/2_trimmed_sanger"
LOG_DIR="$OUT_DIR/logs"
THREADS=8

mkdir -p "$OUT_DIR" "$LOG_DIR"

SAMPLES=$(ls -1 "$IN_DIR"/*-assigned-trustrus_-pair1.fastq 2>/dev/null \
    | sed 's/-assigned-trustrus_-pair1\.fastq//' \
    | xargs -I{} basename {})

echo "Found $(echo "$SAMPLES" | wc -l) assigned sample files. Starting trimming..."

for s in $SAMPLES; do
    echo "  $s"
    skewer \
        -x CCTACGGGNGGCWGCAG \
        -y GGACTACHVGGGTATCTAATCC \
        -m head -k 35 -d 0 -f sanger -t "$THREADS" \
        "$IN_DIR/${s}-assigned-trustrus_-pair1.fastq" \
        "$IN_DIR/${s}-assigned-trustrus_-pair2.fastq" \
        -o "$OUT_DIR/$s" \
        > "$LOG_DIR/${s}.log" 2>&1
done

echo "Compressing output..."
gzip "$OUT_DIR"/*.fastq 2>/dev/null

echo "Done. Logs -> $LOG_DIR/"
