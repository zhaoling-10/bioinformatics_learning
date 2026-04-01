#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bcftools 2>/dev/null || true
cd "$WORKDIR"

mkdir -p variants

echo "[1/2] Calculating SNP statistics..."

# --- Raw VCF stats ---
RAW_TOTAL=$(bcftools view -H variants/all.raw.vcf.gz | wc -l)
RAW_SNPS=$(bcftools view -H -v snps variants/all.raw.vcf.gz | wc -l)
RAW_INDELS=$(bcftools view -H -v indels variants/all.raw.vcf.gz | wc -l)
RAW_SINGLETON=$(bcftools view -H variants/all.raw.vcf.gz \
    | awk '$8 ~ /AC=1;|AC=1$/' | wc -l)

# --- Filtered VCF stats ---
FILT_TOTAL=$(bcftools view -H variants/all.filtered.vcf.gz | wc -l)
FILT_SNPS=$(bcftools view -H -v snps variants/all.filtered.vcf.gz | wc -l)
FILT_INDELS=$(bcftools view -H -v indels variants/all.filtered.vcf.gz | wc -l)
FILT_SINGLETON=$(bcftools view -H variants/all.filtered.vcf.gz \
    | awk '$8 ~ /AC=1;|AC=1$/' | wc -l)

echo "[2/2] Writing summary table..."

cat > variants/snp_summary_figure1.tsv << EOF
category	genome	filter_status	count
total_variants	LE_ref	raw	${RAW_TOTAL}
total_variants	LE_ref	filtered	${FILT_TOTAL}
snps_only	LE_ref	raw	${RAW_SNPS}
snps_only	LE_ref	filtered	${FILT_SNPS}
indels	LE_ref	raw	${RAW_INDELS}
indels	LE_ref	filtered	${FILT_INDELS}
singletons	LE_ref	raw	${RAW_SINGLETON}
singletons	LE_ref	filtered	${FILT_SINGLETON}
EOF

echo ""
echo "====================================================="
echo "  Output: variants/snp_summary_figure1.tsv"
echo ""
echo "  NOTE: After mapping to LT reference genome,"
echo "        add LT rows to this table, then plot"
echo "        Figure 1 comparison bar chart in R."
echo "====================================================="