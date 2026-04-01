#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bedtools 2>/dev/null || true
cd "$WORKDIR"

# mkdir -p variants

# --- Step 1: SAF to BED ---
echo "[1/3] Converting SAF to BED..."

test -f "ref/genes_exon_gene.saf" || {
    echo "ERROR: ref/genes_exon_gene.saf not found"
    exit 1
}

# SAF columns: GeneID(1) Chr(2) Start(3) End(4) Strand(5)
# BED columns: Chr(1) Start_0based(2) End(3) GeneID(4) .(5) Strand(6)
awk 'NR>1 {
    print $2 "\t" ($3-1) "\t" $4 "\t" $1 "\t.\t" $5
}' ref/genes_exon_gene.saf \
| sort -k1,1 -k2,2n \
> ref/genes_exon_gene.bed

echo "    Done: ref/genes_exon_gene.bed"
echo "    Total exon intervals: $(wc -l < ref/genes_exon_gene.bed)"

echo "    First 3 lines of BED:"
head -3 ref/genes_exon_gene.bed | \
    awk '{print "    col1="$1" col2="$2" col3="$3" col4="$4" col5="$5" col6="$6}'
echo ""

# --- Annotation function ---
annotate_vcf() {
    local VCF_IN="$1"
    local TSV_OUT="$2"
    local LABEL="$3"

    echo "  Processing: $LABEL"
    echo "  Input:  $VCF_IN"
    echo "  Output: $TSV_OUT"

    test -f "$VCF_IN" || {
        echo "  ERROR: $VCF_IN not found, skipping"
        return 1
    }

    # Step A: VCF to BED
    # Use only 6 columns to match standard BED
    # col1=CHROM col2=START(0-based) col3=END col4=GeneID(placeholder) col5=QUAL col6=FILTER
    # Store REF and ALT in separate temp file to retrieve later

    # Extract all variant info into a tab-separated temp file
    # Format: CHROM  POS  REF  ALT  QUAL  FILTER
    bcftools view -H "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} {
        print $1, $2, $4, $5, $6, $7
    }' > variants/tmp_vcf_info.txt

    echo "  Total variants: $(wc -l < variants/tmp_vcf_info.txt)"

    # Create BED from VCF (3 columns only for intersection)
    # col1=CHROM col2=START(0-based) col3=END
    awk 'BEGIN{OFS="\t"} {
        print $1, ($2-1), $2
    }' variants/tmp_vcf_info.txt > variants/tmp_vcf_3col.bed

    echo ""
    echo "  VCF BED first 2 lines:"
    head -2 variants/tmp_vcf_3col.bed | \
        awk '{print "    col1="$1" col2="$2" col3="$3}'
    echo ""

    # Step B: bedtools intersect
    # Input A (VCF BED): col1=CHROM col2=START col3=END  (3 cols)
    # Input B (gene BED): col1=CHROM col2=START col3=END col4=GeneID col5=. col6=Strand (6 cols)
    # Output with -wa -wb -loj:
    #   col1-3: from A (VCF BED)
    #   col4-9: from B (gene BED)
    #   So GeneID = col7

    bedtools intersect \
        -a variants/tmp_vcf_3col.bed \
        -b ref/genes_exon_gene.bed \
        -wa -wb -loj \
    > variants/tmp_intersect.txt

    echo "  Intersect output first 2 lines:"
    head -2 variants/tmp_intersect.txt | \
        awk '{
            print "    total_cols="NF
            for(i=1;i<=NF;i++) printf "    col"i"="$i"\n"
        }'
    echo ""

    # Step C: Join intersect result with original VCF info
    # intersect has: col1=CHROM col2=START col3=END col4=CHROM_B col5=START_B col6=END_B col7=GeneID col8=. col9=Strand
    # We need to match back to tmp_vcf_info.txt by CHROM and POS

    # Add line number to both files for joining
    awk '{print NR"\t"$0}' variants/tmp_vcf_info.txt   > variants/tmp_vcf_numbered.txt
    awk '{print NR"\t"$0}' variants/tmp_vcf_3col.bed   > variants/tmp_bed_numbered.txt

    # The intersect output rows correspond 1:1 with input VCF BED rows
    # because -loj keeps all A rows
    # So we can paste them together directly

    paste variants/tmp_intersect.txt variants/tmp_vcf_info.txt \
    | awk 'BEGIN{
        OFS="\t"
        print "GeneID","CHROM","POS","REF","ALT","QUAL","FILTER"
    }
    {
        # From intersect (tmp_intersect.txt):
        # $1=CHROM $2=START $3=END $4=CHROM_B $5=START_B $6=END_B $7=GeneID $8=. $9=Strand
        # From vcf_info (tmp_vcf_info.txt) pasted after:
        # $10=CHROM $11=POS $12=REF $13=ALT $14=QUAL $15=FILTER

        gene  = ($7 == "." ? "intergenic" : $7)
        chrom = $10
        pos   = $11
        ref   = $12
        alt   = $13
        qual  = $14
        filt  = $15

        print gene, chrom, pos, ref, alt, qual, filt
    }' \
    | sort -u \
    > "$TSV_OUT"

    # Cleanup
    rm -f variants/tmp_vcf_info.txt \
          variants/tmp_vcf_3col.bed \
          variants/tmp_vcf_numbered.txt \
          variants/tmp_bed_numbered.txt \
          variants/tmp_intersect.txt

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

    echo "  Output first 5 lines:"
    head -5 "$TSV_OUT"
    echo ""
}

# --- Step 2: Annotate raw VCF ---
echo ""
echo "[2/3] Annotating raw VCF..."
annotate_vcf \
    "variants/all.raw.vcf.gz" \
    "variants/all.raw.annotated_geneid.tsv" \
    "Raw VCF"

# --- Step 3: Annotate filtered VCF ---
echo ""
echo "[3/3] Annotating filtered VCF..."
annotate_vcf \
    "variants/all.filtered.vcf.gz" \
    "variants/all.filtered.annotated_geneid.tsv" \
    "Filtered VCF"

echo ""
echo "====================================================="
echo "  Done!"
echo "  variants/all.raw.annotated_geneid.tsv"
echo "  variants/all.filtered.annotated_geneid.tsv"
echo "====================================================="