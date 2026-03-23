#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load fastp 2>/dev/null || true

echo "Staring run"
echo ""

cd "$WORKDIR"

# -r Do not interpret the backslash \ as an escape character.
while read -r srr; do
  echo "==> fastp $srr"

  # define input and output path
  IN1="fastq/${srr}_1.fastq.gz"
  IN2="fastq/${srr}_2.fastq.gz"
  OUT1="trimmed/${srr}_1.fastq.gz"
  OUT2="trimmed/${srr}_2.fastq.gz"

  # check reference genome
  test -f "$IN1" || { echo "Missing $IN1"; exit 1; }
  test -f "$IN2" || { echo "Missing $IN2"; exit 1; }

  # Skip if already trimmed
  if [[ -f "$OUT1" && -f "$OUT2" ]]; then
    echo "  Trimmed exists, skipping: $srr"
    continue
  fi

  fastp \
    --in1 "$IN1" \
    --in2 "$IN2" \
    --out1 "$OUT1" \
    --out2 "$OUT2" \
    --detect_adapter_for_pe \
    --length_required 50 \
    --thread "$THREADS" \
    --html "qc/${srr}.fastp.html" \
    --json "qc/${srr}.fastp.json" \
    2>&1 | tee -a "logs/05_fastp.log"
    # detect_adapter_for_pe --> Automatically detects and removes adapter sequences in paired-end sequencing.
    # Discard reads shorter than 50.
    # Real-time display and permanent storage of logs

    echo "$srr QC-finish!"
    echo ""
    
done < list.txt

echo "All finished"
echo "fastp done. Reports in $WORKDIR/qc/"