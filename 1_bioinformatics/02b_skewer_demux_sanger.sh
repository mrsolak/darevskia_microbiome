#!/usr/bin/env bash
# Demultiplexing with bare 16S primers (trus/trus_, no inline barcodes) and -f sanger.
# Bare primers are more permissive — inline barcode is left for the trim step to remove.

PROJECT_DIR="/mnt/d/Rwd_claude/darevskia_16S_first"
FILT_DIR="$PROJECT_DIR/raw_files/darevskia_filtered"
ADAPT_PATH="$PROJECT_DIR/raw_files/primers.txt"
MATRIX_PATH="$PROJECT_DIR/raw_files/matrix.txt"
OUT_DIR="$PROJECT_DIR/raw_files/1_demultiplexed_sanger"
LOG_DIR="$OUT_DIR/logs"
THREADS=8

mkdir -p "$OUT_DIR" "$LOG_DIR"

SAMPLES=$(ls -1 "$FILT_DIR" | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//')

echo "Found $(echo "$SAMPLES" | wc -l) samples. Starting demultiplexing with -f sanger..."

for s in $SAMPLES; do
    echo "  $s"
    skewer -x "$ADAPT_PATH" -M "$MATRIX_PATH" -b -m head -k 35 -d 0 -f sanger -t "$THREADS" \
        "$FILT_DIR/${s}_R1.fastq.gz" \
        "$FILT_DIR/${s}_R2.fastq.gz" \
        -o "$OUT_DIR/$s" \
        > "$LOG_DIR/${s}.log" 2>&1
done

echo "Done. Logs -> $LOG_DIR/"
