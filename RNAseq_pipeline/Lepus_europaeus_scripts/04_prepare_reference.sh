#!/usr/bin/env bash
set -euo pipefail
# Load the configuration file, which contains variables such as THREADS, WORKDIR, etc.
source ./config.sh
# Load the STAR aligner and samtools tools
# 2>/dev/null: hide error messages || true: do not interrupt the script even if loading fails (may already be in PATH)
module load samtools 2>/dev/null || true 

# check reference genome
test -f "$WORKDIR/ref_src/$REF_FA_GZ" || { echo "Missing $WORKDIR/ref_src/$REF_FA_GZ"; exit 1; }
test -f "$WORKDIR/ref_src/$REF_GFF_ZIP" || { echo "Missing $WORKDIR/ref_src/$REF_GFF_ZIP"; exit 1; }

# Unzip file
echo "Decompressing reference to $WORKDIR/ref/ ..."
gunzip -c "$WORKDIR/ref_src/$REF_FA_GZ" > "$WORKDIR/ref/genome.fa"
unzip -c "$WORKDIR/ref_src/$REF_GFF_ZIP" > "$WORKDIR/ref/genes.gff"

# creat index (e.g fa.fai)
echo "Indexing FASTA..."
samtools faidx "$WORKDIR/ref/genome.fa"

echo "Reference ready:"
ls -lh "$WORKDIR/ref/genome.fa" "$WORKDIR/ref/genes.gff" "$WORKDIR/ref/genome.fa.fai"
