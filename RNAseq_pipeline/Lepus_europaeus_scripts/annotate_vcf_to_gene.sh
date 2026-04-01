#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bedtools 2>/dev/null || true
cd "$WORKDIR"

mkdir -p variants

# --- Step 1: SAF to BED (0-based coordinates) ---
echo "[1/3] Converting SAF to BED..."

test -f "ref/genes_exon_gene.saf" || {
    echo "ERROR: ref/genes_exon_gene.saf not found"
    exit 1
}

awk 'NR>1 {
    print $2 "\t" ($3-1) "\t" $4 "\t" $1 "\t.\t" $5
}' ref/genes_exon_gene.saf \
| sort -k1,1 -k2,2n \
> ref/genes_exon_gene.bed

echo "    Done: ref/genes_exon_gene.bed"
echo "    Total exon intervals: $(wc -l < ref/genes_exon_gene.bed)"

# --- Annotation function ---
annotate_vcf() {
    local VCF_IN="$1"
    local TSV_OUT="$2"
    local LABEL="$3"

    echo ""
    echo "  Processing: $LABEL"
    echo "  Input:  $VCF_IN"
    echo "  Output: $TSV_OUT"

    test -f "$VCF_IN" || {
        echo "  ERROR: $VCF_IN not found, skipping"
        return 1
    }

    bcftools view -H "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} {
        chrom=$1; pos=$2; ref=$4; alt=$5; qual=$6; filt=$7
        print chrom, (pos-1), pos, chrom"_"pos"_"ref"_"alt, qual, filt
    }' > variants/tmp_vcf.bed

    echo "  Total variant sites: $(wc -l < variants/tmp_vcf.bed)"

    bedtools intersect \
        -a variants/tmp_vcf.bed \
        -b ref/genes_exon_gene.bed \
        -wa -wb -loj \
    | awk 'BEGIN{
        OFS="\t"
        print "GeneID","CHROM","POS","REF","ALT","QUAL","FILTER"
    }
    NR>0 {
        chrom = $1
        pos   = $3
        qual  = $5
        filt  = $6
        name  = $4
        n = split(name, parts, "_")
        alt = parts[n]
        ref = parts[n-1]
        gene = ($11 == "." ? "intergenic" : $11)
        print gene, chrom, pos, ref, alt, qual, filt
    }' \
    | sort -u \
    > "$TSV_OUT"

    rm -f variants/tmp_vcf.bed

    local TOTAL IN_GENE INTERGENIC
    TOTAL=$(tail -n +2 "$TSV_OUT" | wc -l)
    IN_GENE=$(tail -n +2 "$TSV_OUT" | awk '$1 != "intergenic"' | wc -l)
    INTERGENIC=$(tail -n +2 "$TSV_OUT" | awk '$1 == "intergenic"' | wc -l)

    echo "  --- Stats ---"
    echo "  Total annotated rows:    $TOTAL"
    echo "  Variants in gene regions: $IN_GENE"
    echo "  Intergenic variants:      $INTERGENIC"
    if [ "$TOTAL" -gt 0 ]; then
        awk -v ig="$IN_GENE" -v t="$TOTAL" \
            'BEGIN{printf "  Gene region coverage:     %.1f%%\n", ig/t*100}'
    fi
    echo "  -------------"
}

# --- Step 2: Annotate raw VCF ---
echo ""
echo "[2/3] Annotating raw VCF (all.raw.vcf.gz)..."
annotate_vcf \
    "variants/all.raw.vcf.gz" \
    "variants/all.raw.annotated_geneid.tsv" \
    "Raw VCF"

# --- Step 3: Annotate filtered VCF ---
echo ""
echo "[3/3] Annotating filtered VCF (all.filtered.vcf.gz)..."
annotate_vcf \
    "variants/all.filtered.vcf.gz" \
    "variants/all.filtered.annotated_geneid.tsv" \
    "Filtered VCF"

echo ""
echo "====================================================="
echo "  Output files:"
echo "  variants/all.raw.annotated_geneid.tsv"
echo "  variants/all.filtered.annotated_geneid.tsv"
echo "====================================================="