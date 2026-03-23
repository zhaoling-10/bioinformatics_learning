#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load samtools 2>/dev/null || true

echo "start Post-alignment QC"
echo ""

cd "$WORKDIR"

# flagstat (The following section is running flagstat quality control)
while read -r srr; do
  BAM="aln/${srr}.Aligned.sortedByCoord.out.bam"
  test -f "$BAM" || { echo "Missing $BAM"; exit 1; }
  samtools flagstat "$BAM" > "qc/${srr}.flagstat.txt"
done < list.txt
# Define BAM file path
# Check if the BAM file exists
# Core command: run samtools flagstat(Function of flagstat: summarize alignment information in a BAM file)

# STAR mapping summary
# echo -e: enable escape character interpretation, create a new file in qc and write the header "sampletuniquely_mapped_pct"
echo -e "sample\tuniquely_mapped_pct" > qc/star_mapping_summary.tsv
while read -r srr; do
  LOG="aln/${srr}.Log.final.out"
  # define file path
  test -f "$LOG" || { echo "Missing $LOG"; exit 1; }
  # check whether the files exists
  # Extract the unique match rate from the log, remove the spaces before and after, and copy it to 'uniq'
  uniq=$(grep "Uniquely mapped reads %" "$LOG" | awk -F'|\t' '{print $2}' | xargs)
  # Append the results to the summary file without overwriting existing content, adding a new line at the end of the file
  echo -e "${srr}\t${uniq}" >> qc/star_mapping_summary.tsv
done < list.txt

echo "Post-alignment QC done."