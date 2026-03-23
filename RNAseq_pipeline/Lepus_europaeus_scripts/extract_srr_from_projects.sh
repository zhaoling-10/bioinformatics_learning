#!/bin/bash
set -euo pipefail

BASE_DIR="/scratch/project_2002674/RNAseq_hares"
METADATA_DIR="$BASE_DIR/metadata"
RUNINFO_DIR="$METADATA_DIR/runinfo"
PROJECT_LIST="$METADATA_DIR/bioproject_list.txt"
ALL_SRR_FILE="$METADATA_DIR/all_SRR.txt"
MAPPING_FILE="$METADATA_DIR/project_to_srr_mapping.txt"  # ★ Add mapping file

mkdir -p "$RUNINFO_DIR"

# Clear old files
> "$ALL_SRR_FILE"
> "$MAPPING_FILE"  # ★ Clear mapping file

# Add mapping file header
echo -e "BioProject\tSRR\tStudy_Title\tSample_Name" > "$MAPPING_FILE"

echo "============================================"
echo "Extracting SRR from BioProjects"
echo "Start: $(date)"
echo "============================================"
echo ""

while read PROJECT; do
    [[ -z "$PROJECT" ]] && continue
    [[ "$PROJECT" =~ ^# ]] && continue
    
    echo "=============================="
    echo "Processing $PROJECT"
    date
    
    OUT_TSV="$RUNINFO_DIR/${PROJECT}_ena_runs.tsv"
    
    # Retrieve data from the ENA API (containing more information)
    curl -s \
      "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession=${PROJECT}&fields=run_accession,study_title,sample_title,sample_alias&format=tsv&limit=0" \
      -o "$OUT_TSV"
    
    # Extract SRR and add it to the overall list
    awk 'NR>1 && $1 ~ /^[SDE]RR/ {print $1}' "$OUT_TSV" >> "$ALL_SRR_FILE"
    
    # ★★★ Create mapping: Project ID -> SRR -> Sample Information ★★★
    awk -v proj="$PROJECT" 'NR>1 && $1 ~ /^[SDE]RR/ {print proj"\t"$1"\t"$2"\t"$3}' "$OUT_TSV" >> "$MAPPING_FILE"
    
    echo "Finished $PROJECT"
    echo "SRR count: $(awk 'NR>1 && $1 ~ /^[SDE]RR/' "$OUT_TSV" | wc -l)"
    date
    echo ""
    
done < "$PROJECT_LIST"

# Deduplication SRR list
sort -u "$ALL_SRR_FILE" -o "$ALL_SRR_FILE"

echo ""
echo "=============================="
echo "SUMMARY"
echo "=============================="
echo "Total unique SRRs: $(wc -l < $ALL_SRR_FILE)"
echo ""
echo "Files created:"
echo "  SRR list: $ALL_SRR_FILE"
echo "  Mapping file: $MAPPING_FILE"
echo ""
echo "Mapping file preview:"
head -10 "$MAPPING_FILE"
echo ""
echo "Finished: $(date)"
echo "=============================="
