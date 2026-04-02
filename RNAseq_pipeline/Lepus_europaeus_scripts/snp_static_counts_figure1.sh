#!/usr/bin/env bash
set -euo pipefail

module load bcftools 2>/dev/null || true
module load bedtools 2>/dev/null || true

LE_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/lepus_europaeus_gff"
LT_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/lepus_timidus_gff"
OUT_DIR="${LE_DIR}/figure1_data"

mkdir -p "$OUT_DIR"

# --- Step 1: Check all input files ---
echo "[0/5] Checking input files..."

ALL_OK=true
for f in \
    "${LE_DIR}/variants/all.filtered.vcf.gz" \
    "${LT_DIR}/variants/all.filtered.vcf.gz" \
    "${LE_DIR}/ref/genes_exon_gene.saf" \
    "${LT_DIR}/ref/genes_exon_gene.saf"; do
    if [ -f "$f" ]; then
        echo "  FOUND:   $f"
    else
        echo "  MISSING: $f"
        ALL_OK=false
    fi
done

# Check raw VCF separately (optional for Figure 1C)
for f in \
    "${LE_DIR}/variants/all.raw.vcf.gz" \
    "${LT_DIR}/variants/all.raw.vcf.gz"; do
    if [ -f "$f" ]; then
        echo "  FOUND:   $f"
    else
        echo "  WARNING: $f not found (needed for raw SNP counts)"
    fi
done

$ALL_OK || {
    echo ""
    echo "ERROR: Missing required files. Please check paths."
    exit 1
}
echo ""

# --- Step 2: SAF to BED ---
# SAF format: GeneID(1) Chr(2) Start_1based(3) End(4) Strand(5)
# BED format: Chr(1) Start_0based(2) End(3) GeneID(4) .(5) Strand(6)
# Key: BED_Start = SAF_Start - 1  (coordinate system conversion)
# Key: BED has no header line (bedtools requirement)
# Key: BED must be sorted by Chr then Start (bedtools -sorted requirement)

echo "[1/5] Converting SAF to BED..."

awk 'NR>1 {
    print $2 "\t" ($3-1) "\t" $4 "\t" $1 "\t.\t" $5
}' "${LE_DIR}/ref/genes_exon_gene.saf" \
| sort -k1,1 -k2,2n \
> "${LE_DIR}/ref/genes_exon_gene.bed"

awk 'NR>1 {
    print $2 "\t" ($3-1) "\t" $4 "\t" $1 "\t.\t" $5
}' "${LT_DIR}/ref/genes_exon_gene.saf" \
| sort -k1,1 -k2,2n \
> "${LT_DIR}/ref/genes_exon_gene.bed"

echo "  LE BED: $(wc -l < ${LE_DIR}/ref/genes_exon_gene.bed) intervals"
echo "  LT BED: $(wc -l < ${LT_DIR}/ref/genes_exon_gene.bed) intervals"

# Verify BED conversion
echo ""
echo "  SAF vs BED comparison (LE, first gene):"
echo "  SAF: $(grep -m1 'NC_084827.1-g1' ${LE_DIR}/ref/genes_exon_gene.saf)"
echo "  BED: $(grep -m1 'NC_084827.1-g1' ${LE_DIR}/ref/genes_exon_gene.bed)"
echo ""

# --- Step 3: SNP statistics for Figure 1C ---
echo "[2/5] Calculating SNP statistics for Figure 1C..."

get_snp_stats() {
    local VCF="$1"
    local GENOME="$2"
    local FILTER="$3"

    test -f "$VCF" || {
        echo "  WARNING: $VCF not found, skipping $GENOME $FILTER"
        return
    }

    local SNPS SINGLETONS
    SNPS=$(bcftools view -H -v snps "$VCF" | wc -l)
    SINGLETONS=$(bcftools view -H -v snps "$VCF" | \
        awk '$8 ~ /AC=1(;|$)/' | wc -l)

    echo -e "snps\t${GENOME}\t${FILTER}\t${SNPS}"
    echo -e "singletons\t${GENOME}\t${FILTER}\t${SINGLETONS}"
}

{
    echo -e "category\tgenome\tfilter_status\tcount"
    get_snp_stats \
        "${LE_DIR}/variants/all.raw.vcf.gz"      "LE_ref" "raw"
    get_snp_stats \
        "${LE_DIR}/variants/all.filtered.vcf.gz" "LE_ref" "filtered"
    get_snp_stats \
        "${LT_DIR}/variants/all.raw.vcf.gz"      "LT_ref" "raw"
    get_snp_stats \
        "${LT_DIR}/variants/all.filtered.vcf.gz" "LT_ref" "filtered"
} > "${OUT_DIR}/snp_summary.tsv"

echo "  Done: ${OUT_DIR}/snp_summary.tsv"
echo "  Content:"
cat "${OUT_DIR}/snp_summary.tsv"
echo ""

# --- Step 4: Annotate VCFs with gene IDs ---
echo "[3/5] Annotating VCFs with gene IDs..."

annotate_vcf() {
    local VCF_IN="$1"
    local BED="$2"
    local TSV_OUT="$3"
    local LABEL="$4"
    local TMP="${OUT_DIR}/tmp_${LABEL}"

    echo "  Processing: $LABEL"
    mkdir -p "$TMP"

    # Extract: CHROM POS REF ALT QUAL FILTER
    bcftools view -H "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} {
        print $1, $2, $4, $5, $6, $7
    }' > "${TMP}/vcf_info.txt"

    echo "  Total variants: $(wc -l < ${TMP}/vcf_info.txt)"

    # VCF to 3-col sorted BED
    awk 'BEGIN{OFS="\t"} {
        print $1, ($2-1), $2
    }' "${TMP}/vcf_info.txt" \
    | sort -k1,1 -k2,2n \
    > "${TMP}/vcf.bed"

    # Gene BED already sorted, just copy
    sort -k1,1 -k2,2n "$BED" > "${TMP}/gene.bed"

    # Variants IN gene regions
    bedtools intersect \
        -a "${TMP}/vcf.bed" \
        -b "${TMP}/gene.bed" \
        -wa -wb \
        -sorted \
    > "${TMP}/in_gene.bed"

    # Intergenic variants
    bedtools intersect \
        -a "${TMP}/vcf.bed" \
        -b "${TMP}/gene.bed" \
        -v \
        -sorted \
    > "${TMP}/intergenic.bed"

    echo "  In-gene variants:   $(wc -l < ${TMP}/in_gene.bed)"
    echo "  Intergenic variants: $(wc -l < ${TMP}/intergenic.bed)"

    # Build TSV output
    echo -e "GeneID\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER" \
        > "$TSV_OUT"

    # in_gene.bed columns:
    # $1=CHROM $2=START $3=END(=POS) $4=CHROM_B $5=START_B $6=END_B $7=GeneID $8=. $9=Strand
    if [ -s "${TMP}/in_gene.bed" ]; then
        awk 'BEGIN{OFS="\t"}
        NR==FNR {
            key = $1"_"$2
            info[key] = $3"\t"$4"\t"$5"\t"$6
            next
        }
        {
            key = $1"_"$3
            gene = $7
            if (key in info)
                print gene, $1, $3, info[key]
        }' "${TMP}/vcf_info.txt" \
           "${TMP}/in_gene.bed" \
        | sort -u >> "$TSV_OUT"
    fi

    if [ -s "${TMP}/intergenic.bed" ]; then
        awk 'BEGIN{OFS="\t"}
        NR==FNR {
            key = $1"_"$2
            info[key] = $3"\t"$4"\t"$5"\t"$6
            next
        }
        {
            key = $1"_"$3
            if (key in info)
                print "intergenic", $1, $3, info[key]
        }' "${TMP}/vcf_info.txt" \
           "${TMP}/intergenic.bed" \
        | sort -u >> "$TSV_OUT"
    fi

    rm -rf "$TMP"

    local IN_GENE INTERGENIC TOTAL
    TOTAL=$(tail -n +2 "$TSV_OUT" | wc -l)
    IN_GENE=$(awk 'NR>1 && $1!="intergenic"' "$TSV_OUT" | wc -l)
    INTERGENIC=$(awk 'NR>1 && $1=="intergenic"' "$TSV_OUT" | wc -l)

    echo "  Output rows total:  $TOTAL"
    echo "  Gene variants:      $IN_GENE"
    echo "  Intergenic:         $INTERGENIC"
    if [ "$TOTAL" -gt 0 ]; then
        awk -v ig="$IN_GENE" -v t="$TOTAL" \
            'BEGIN{printf "  Gene coverage:      %.1f%%\n", ig/t*100}'
    fi
    echo "  First 3 output lines:"
    head -4 "$TSV_OUT" | column -t
    echo ""
}

annotate_vcf \
    "${LE_DIR}/variants/all.filtered.vcf.gz" \
    "${LE_DIR}/ref/genes_exon_gene.bed" \
    "${OUT_DIR}/LE_annotated.tsv" \
    "LE"

annotate_vcf \
    "${LT_DIR}/variants/all.filtered.vcf.gz" \
    "${LT_DIR}/ref/genes_exon_gene.bed" \
    "${OUT_DIR}/LT_annotated.tsv" \
    "LT"

# --- Step 5: Gene-level SNP overlap for Figure 1D ---
echo "[4/5] Calculating gene-level SNP overlap for Figure 1D..."

awk 'NR>1 && $1!="intergenic" {print $1}' \
    "${OUT_DIR}/LE_annotated.tsv" \
| sort -u > "${OUT_DIR}/LE_genes.txt"

awk 'NR>1 && $1!="intergenic" {print $1}' \
    "${OUT_DIR}/LT_annotated.tsv" \
| sort -u > "${OUT_DIR}/LT_genes.txt"

comm -12 "${OUT_DIR}/LE_genes.txt" \
         "${OUT_DIR}/LT_genes.txt" \
> "${OUT_DIR}/shared_genes.txt"

comm -23 "${OUT_DIR}/LE_genes.txt" \
         "${OUT_DIR}/LT_genes.txt" \
> "${OUT_DIR}/LE_only_genes.txt"

comm -13 "${OUT_DIR}/LE_genes.txt" \
         "${OUT_DIR}/LT_genes.txt" \
> "${OUT_DIR}/LT_only_genes.txt"

LE_TOTAL=$(wc -l < "${OUT_DIR}/LE_genes.txt")
LT_TOTAL=$(wc -l < "${OUT_DIR}/LT_genes.txt")
SHARED=$(wc -l   < "${OUT_DIR}/shared_genes.txt")
LE_ONLY=$(wc -l  < "${OUT_DIR}/LE_only_genes.txt")
LT_ONLY=$(wc -l  < "${OUT_DIR}/LT_only_genes.txt")

cat > "${OUT_DIR}/snp_overlap_summary.tsv" << EOF
snp_type	genome	count
shared	LE_ref	${SHARED}
shared	LT_ref	${SHARED}
unique	LE_ref	${LE_ONLY}
unique	LT_ref	${LT_ONLY}
EOF

echo "  LE genes with SNPs: $LE_TOTAL"
echo "  LT genes with SNPs: $LT_TOTAL"
echo "  Shared:             $SHARED"
echo "  LE only:            $LE_ONLY"
echo "  LT only:            $LT_ONLY"
echo ""

echo "[5/5] All done!"
echo ""
echo "============================================"
echo "  Output directory: $OUT_DIR"
echo "  snp_summary.tsv          -> Figure 1C"
echo "  snp_overlap_summary.tsv  -> Figure 1D"
echo "  LE_annotated.tsv         -> LE variants"
echo "  LT_annotated.tsv         -> LT variants"
echo "============================================"
echo "  Next: run Rscript plot_figure1CD.R"
echo "============================================"