#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load subread 2>/dev/null || true

cd "$WORKDIR"

echo "start gene_counts_featureCounts"

# check whether GTF files exist
test -f "ref/genes.gtf" || { echo "Missing"; exit 1; }

# Use the featureCounts command (from the subread package)
# This is currently the most commonly used gene counting tool, fast and with low memory usage
featureCounts \
  -T "$THREADS" \
  -p \
  -s 0 \
  -a ref/genes.gtf \
  -o counts/gene_counts_featureCounts.txt \
  aln/*.Aligned.sortedByCoord.out.bam \
  2>&1 | tee logs/09_featurecounts.log

echo "featureCounts done: counts/gene_counts_featureCounts.txt"
echo "NOTE: -s 0 assumes unstranded library. If stranded, use -s 1 or -s 2."
# indicates that the input is paired-end data; Tells featureCounts to count a pair of reads as a single fragment
# Set library strand specificity
# -s 0: Non-strand-specific library (unstranded)
# -s 1: Strand-specific, read1 is in the same direction as the gene
# -s 2: Strand-specific, read1 is in the opposite direction to the gene
# Specified gene annotation file 
# Specify output file
# Input File (output file and show in screen)