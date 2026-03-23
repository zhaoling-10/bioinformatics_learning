#!/bin/bash
#SBATCH --job-name=ENA_download
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2002674/RNAseq_hares/logs/ena_%j.out
#SBATCH --error=/scratch/project_2002674/RNAseq_hares/logs/ena_%j.err

echo "============================================"
echo "Download RNA-seq data using enaDataGet"
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "============================================"
echo ""

module load biokit

WORK_DIR="/scratch/project_2002674/RNAseq_hares"
SRR_LIST="${WORK_DIR}/metadata/all_SRR.txt"
FASTQ_DIR="${WORK_DIR}/fastq"

mkdir -p "$FASTQ_DIR"

# Deduplication SRR list
sort -u "$SRR_LIST" > "${WORK_DIR}/metadata/all_SRR_unique.txt"
SRR_LIST="${WORK_DIR}/metadata/all_SRR_unique.txt"

TOTAL=$(wc -l < "$SRR_LIST")
echo "Total unique samples: $TOTAL"
echo ""

CURRENT=0
SUCCESS=0
FAILED=0
SKIPPED=0

while IFS= read -r SRR || [[ -n "$SRR" ]]; do
    [[ -z "$SRR" ]] && continue
    [[ "$SRR" =~ ^# ]] && continue
    
    CURRENT=$((CURRENT + 1))
    
    echo "=========================================="
    echo "[$CURRENT/$TOTAL] $SRR"
    echo "Time: $(date '+%H:%M:%S')"
    echo "=========================================="
    
    # ★★★ Key: Check if it has been downloaded ★★★
    if [[ -d "$FASTQ_DIR/$SRR" ]]; then
        # Check if the folder contains content.
        FILE_COUNT=$(ls "$FASTQ_DIR/$SRR"/*.fastq.gz 2>/dev/null | wc -l)
        if [[ $FILE_COUNT -gt 0 ]]; then
            echo "  ⊙ Already downloaded, skipping..."
            SKIPPED=$((SKIPPED + 1))
            echo ""
            continue
        fi
    fi
    
    # download dataset
    if enaDataGet -f fastq -d "$FASTQ_DIR" "$SRR" 2>&1; then
        echo "  ✓ Downloaded successfully"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "  ✗ Download FAILED"
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
    
done < "$SRR_LIST"

echo ""
echo "============================================"
echo "DOWNLOAD SUMMARY"
echo "============================================"
echo "Total:    $CURRENT"
echo "New downloads: $SUCCESS"
echo "Skipped (already exists): $SKIPPED"
echo "Failed:   $FAILED"
echo ""

TOTAL_FILES=$(ls -d "$FASTQ_DIR"/SRR*/ 2>/dev/null | wc -l)
echo "Total samples in directory: $TOTAL_FILES"
echo "Total size: $(du -sh $FASTQ_DIR | cut -f1)"
echo ""

echo "Finished: $(date)"
echo "============================================"

exit 0
