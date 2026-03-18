#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load star 2>/dev/null || true
module load samtools 2>/dev/null || true

# check reference genome
test -f "$WORKDIR/ref/genome.fa" || { echo "Missing ref/genome.fa (run 04_prepare_reference.sh)"; exit 1; }
test -f "$WORKDIR/ref/genes.gff" || { echo "Missing ref/genes.gff (run 04_prepare_reference.sh)"; exit 1; }

cd "$WORKDIR"

# STAR index is created once
if [ -f "$WORKDIR/star_index/Genome" ]; then
  echo "STAR index already exists. Skipping."
  exit 0
fi

STAR \
  --runMode genomeGenerate \                        # Mode: generate genome index
  --runThreadN "$THREADS" \                         # Number of threads to use
  --genomeDir "$WORKDIR/star_index" \               # Index output directory
  --genomeFastaFiles "$WORKDIR/ref/genome.fa" \     # Reference genome
  --sjdbGTFfile "$WORKDIR/ref/genes.gff" \          # Gene annotation (including splice sites)
  --sjdbOverhang "$SJDB_OVERHANG" \                 # Read length minus 1 (e.g., 100)
  2>&1 | tee "$WORKDIR/logs/06_star_index.log"

echo "STAR index done."
# need big capacity and long time