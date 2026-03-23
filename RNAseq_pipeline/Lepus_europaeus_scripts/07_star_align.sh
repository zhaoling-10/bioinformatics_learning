#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load star 2>/dev/null || true
module load samtools 2>/dev/null || true

echo "start STAR align..."
echo ""

cd "$WORKDIR"

# Loop: Read SRR numbers line by line from the list.txt file
while read -r srr; do
  echo "==> STAR align $srr"

  # definition and check files(PATH AND check)
  R1="trimmed/${srr}_1.fastq.gz"
  R2="trimmed/${srr}_2.fastq.gz"
  test -f "$R1" || { echo "Missing $R1"; exit 1; }
  test -f "$R2" || { echo "Missing $R2"; exit 1; }

  # definition outcome file(The file name format is the default of STAR, indicating a BAM file sorted by coordinates)
  BAM="aln/${srr}.Aligned.sortedByCoord.out.bam"
  # Only run the alignment if the BAM file does not exist
  if [ ! -f "$BAM" ]; then
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$WORKDIR/star_index" \
      --readFilesIn "$R1" "$R2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "aln/${srr}." \
      --outSAMtype BAM SortedByCoordinate \
      2>&1 | tee -a "logs/07_star_align.log"
  else
    echo "  BAM exists, skipping: $srr"
  fi

  # create BAM index (Check if the BAM index file exists, and if it does not, create the index (.bai file) using samtools index)
  if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
  fi
done < list.txt

echo "STAR alignment done."
# runThreadN "$THREADS": umber of threads to use (read from config.sh)
# genomeDir "$WORKDIR/star_index": Specify the previously built index 
# readFilesIn "$R1" "$R2": Input paired-end FASTQ files
# readFilesCommand zcat: Command to read compressed files
# outFileNamePrefix "aln/${srr}.": Prefix for output files (including path)
# outSAMtype BAM SortedByCoordinate: SortedByCoordinate Output BAM file sorted by coordinates
# 2>&1 | tee -a "logs/07_star_align.log": Append both standard error and standard output to the log file while displaying on screen