#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bedtools 2>/dev/null || true
cd "$WORKDIR"

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

# --- DEBUG: check actual column structure ---
echo ""
echo "    DEBUG - BED file first 3 lines:"
head -3 ref/genes_exon_gene.bed
echo ""

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

    # VCF to BED
    # Format: CHROM  START(0-based)  END  CHROM  POS  REF  ALT  QUAL  FILTER
    # Keep REF and ALT as separate fields to avoid underscore parsing issues
    bcftools view -H "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} {
        chrom=$1; pos=$2; ref=$4; alt=$5; qual=$6; filt=$7
        print chrom, (pos-1), pos, chrom, pos, ref, alt, qual, filt
    }' > variants/tmp_vcf.bed

    echo "  Total variant sites: $(wc -l < variants/tmp_vcf.bed)"

    # DEBUG: check VCF BED structure
    echo "  DEBUG - VCF BED first 2 lines:"
    head -2 variants/tmp_vcf.bed
    echo ""

    # Run bedtools intersect
    # Input VCF BED columns:  1=CHROM 2=START 3=END 4=CHROM 5=POS 6=REF 7=ALT 8=QUAL 9=FILTER
    # Input gene BED columns: 1=CHROM 2=START 3=END 4=GeneID 5=. 6=Strand
    # After -wa -wb -loj:
    # Cols 1-9:  from VCF BED
    # Cols 10-15: from gene BED
    # So GeneID = col 13

    bedtools intersect \
        -a variants/tmp_vcf.bed \
        -b ref/genes_exon_gene.bed \
        -wa -wb -loj \
    > variants/tmp_intersect.bed

    # DEBUG: check intersect output structure
    echo "  DEBUG - Intersect output first 2 lines:"
    head -2 variants/tmp_intersect.bed
    echo "  DEBUG - Total columns in intersect output:"
    head -1 variants/tmp_intersect.bed | awk '{print NF}'
    echo ""

    # Parse intersect output
    # VCF BED:  col1=CHROM col2=START col3=END col4=CHROM col5=POS col6=REF col7=ALT col8=QUAL col9=FILTER
    # gene BED: col10=CHROM col11=START col12=END col13=GeneID col14=. col15=Strand
    awk 'BEGIN{
        OFS="\t"
        print "GeneID","CHROM","POS","REF","ALT","QUAL","FILTER"
    }
    {
        chrom = $1
        pos   = $5
        ref   = $6
        alt   = $7
        qual  = $8
        filt  = $9
        gene  = ($13 == "." ? "intergenic" : $13)
        print gene, chrom, pos, ref, alt, qual, filt
    }' variants/tmp_intersect.bed \
    | sort -u \
    > "$TSV_OUT"

#    rm -f variants/tmp_vcf.bed variants/tmp_intersect.bed

    # Stats
    local TOTAL IN_GENE INTERGENIC
    TOTAL=$(tail -n +2 "$TSV_OUT" | wc -l)
    IN_GENE=$(tail -n +2 "$TSV_OUT" | awk '$1 != "intergenic"' | wc -l)
    INTERGENIC=$(tail -n +2 "$TSV_OUT" | awk '$1 == "intergenic"' | wc -l)

    echo "  --- Stats ---"
    echo "  Total annotated rows:     $TOTAL"
    echo "  Variants in gene regions: $IN_GENE"
    echo "  Intergenic variants:      $INTERGENIC"
    if [ "$TOTAL" -gt 0 ]; then
        awk -v ig="$IN_GENE" -v t="$TOTAL" \
            'BEGIN{printf "  Gene region coverage:     %.1f%%\n", ig/t*100}'
    fi
    echo "  -------------"

    # DEBUG: show first 5 lines of output
    echo "  DEBUG - Output first 5 lines:"
    head -5 "$TSV_OUT"
    echo ""
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