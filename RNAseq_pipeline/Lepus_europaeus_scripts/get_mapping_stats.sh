#!/usr/bin/env bash
set -euo pipefail

LE_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/lepus_europaeus_gff"
LT_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/lepus_timidus_gff"
OUT="${LE_DIR}/figure1_data/mapping_stats.tsv"

mkdir -p "${LE_DIR}/figure1_data"

echo -e "sample\tgenome\tspecies\tuniquely_mapped_pct\ttotal_reads\tmapped_reads" > "$OUT"

# LE species samples
LE_SAMPLES="SRR18740835 SRR18740836 SRR18740837 SRR18740838"
# LT species samples
LT_SAMPLES="SRR18740839 SRR18740840 SRR18740841 SRR18740842"

parse_log() {
    local LOG="$1"
    local SAMPLE="$2"
    local GENOME="$3"
    local SPECIES="$4"

    local UNIQ_PCT TOTAL MAPPED

    UNIQ_PCT=$(grep "Uniquely mapped reads %" "$LOG" \
        | awk -F'|' '{gsub(/[[:space:]]|%/,"",$2); print $2}')

    TOTAL=$(grep "Number of input reads" "$LOG" \
        | awk -F'|' '{gsub(/[[:space:]]/,"",$2); print $2}')

    MAPPED=$(grep "Uniquely mapped reads number" "$LOG" \
        | awk -F'|' '{gsub(/[[:space:]]/,"",$2); print $2}')

    echo -e "${SAMPLE}\t${GENOME}\t${SPECIES}\t${UNIQ_PCT}\t${TOTAL}\t${MAPPED}"
}

echo "Processing LE samples mapped to LE reference..."
for SAMPLE in $LE_SAMPLES; do
    LOG="${LE_DIR}/aln/${SAMPLE}.Log.final.out"
    parse_log "$LOG" "$SAMPLE" "LE_ref" "LE"
done >> "$OUT"

echo "Processing LT samples mapped to LE reference..."
for SAMPLE in $LT_SAMPLES; do
    LOG="${LE_DIR}/aln/${SAMPLE}.Log.final.out"
    parse_log "$LOG" "$SAMPLE" "LE_ref" "LT"
done >> "$OUT"

echo "Processing LE samples mapped to LT reference..."
for SAMPLE in $LE_SAMPLES; do
    LOG="${LT_DIR}/aln/${SAMPLE}.Log.final.out"
    parse_log "$LOG" "$SAMPLE" "LT_ref" "LE"
done >> "$OUT"

echo "Processing LT samples mapped to LT reference..."
for SAMPLE in $LT_SAMPLES; do
    LOG="${LT_DIR}/aln/${SAMPLE}.Log.final.out"
    parse_log "$LOG" "$SAMPLE" "LT_ref" "LT"
done >> "$OUT"

echo ""
echo "Done: $OUT"
echo ""
echo "Content:"
cat "$OUT"