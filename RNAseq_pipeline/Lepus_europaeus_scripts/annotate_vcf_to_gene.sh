#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load biokit 2>/dev/null || \
module load bedtools 2>/dev/null || \
module load bcftools 2>/dev/null || true

cd "$WORKDIR"

# mkdir -p variants

# --- Step 1: Check tools ---
echo "[0/4] Checking required tools..."

for tool in bcftools bedtools; do
    if ! command -v $tool &>/dev/null; then
        echo "ERROR: $tool not found. Trying to load..."
        module load $tool 2>/dev/null || {
            echo "ERROR: Cannot load $tool. Please run:"
            echo "  module load bcftools"
            echo "  module load bedtools"
            exit 1
        }
    fi
    echo "  $tool: $(command -v $tool)"
done
echo ""

# --- Step 2: SAF to BED ---
echo "[1/4] Converting SAF to BED..."

test -f "ref/genes_exon_gene.saf" || {
    echo "ERROR: ref/genes_exon_gene.saf not found"
    exit 1
}

# SAF: GeneID(1) Chr(2) Start_1based(3) End(4) Strand(5)
# BED: Chr(1) Start_0based(2) End(3) GeneID(4) .(5) Strand(6)
awk 'NR>1 {
    print $2 "\t" ($3-1) "\t" $4 "\t" $1 "\t.\t" $5
}' ref/genes_exon_gene.saf \
| sort -k1,1 -k2,2n \
> ref/genes_exon_gene.bed

echo "    Done: ref/genes_exon_gene.bed"
echo "    Total exon intervals: $(wc -l < ref/genes_exon_gene.bed)"
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

    # Step A: Extract VCF info to temp file
    # Format: CHROM  POS  REF  ALT  QUAL  FILTER
    bcftools view -H "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} {
        print $1, $2, $4, $5, $6, $7
    }' > variants/tmp_vcf_info.txt

    local NVAR
    NVAR=$(wc -l < variants/tmp_vcf_info.txt)
    echo "  Total variants: $NVAR"

    # Step B: Convert VCF to 3-column BED for intersection
    # CHROM  START(0-based)  END
    awk 'BEGIN{OFS="\t"} {
        print $1, ($2-1), $2
    }' variants/tmp_vcf_info.txt \
    | sort -k1,1 -k2,2n \
    > variants/tmp_vcf_sorted.bed

    # Step C: Sort gene BED (already done but re-sort to be safe)
    sort -k1,1 -k2,2n ref/genes_exon_gene.bed \
    > variants/tmp_gene_sorted.bed

    # Step D: Use bedtools intersect WITHOUT -loj
    # Run two separate intersections:
    # 1. Find variants that ARE in gene regions
    # 2. Subtract to find variants NOT in gene regions

    # Intersection: variants IN gene regions
    # Output: all columns from A + all columns from B
    # A cols (3): CHROM START END
    # B cols (6): CHROM START END GeneID . Strand
    # Total 9 cols, GeneID = col7
    bedtools intersect \
        -a variants/tmp_vcf_sorted.bed \
        -b variants/tmp_gene_sorted.bed \
        -wa -wb \
        -sorted \
    > variants/tmp_in_gene.bed

    echo "  Variants overlapping gene regions: $(wc -l < variants/tmp_in_gene.bed)"

    # Verify column structure
    if [ -s variants/tmp_in_gene.bed ]; then
        echo "  Column check (first line of gene-overlapping variants):"
        head -1 variants/tmp_in_gene.bed | \
            awk '{for(i=1;i<=NF;i++) printf "    col%d=%s\n",i,$i}'
        echo ""
    fi

    # Variants NOT in gene regions
    bedtools intersect \
        -a variants/tmp_vcf_sorted.bed \
        -b variants/tmp_gene_sorted.bed \
        -v \
        -sorted \
    > variants/tmp_intergenic.bed

    echo "  Intergenic variants: $(wc -l < variants/tmp_intergenic.bed)"
    echo ""

    # Step E: Build output TSV

    # Write header
    echo -e "GeneID\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER" > "$TSV_OUT"

    # Process in-gene variants
    # col1=CHROM col2=START col3=END col4=CHROM_B col5=START_B col6=END_B col7=GeneID col8=. col9=Strand
    # Match back to vcf_info by CHROM and POS (POS = END in 3-col BED)
    if [ -s variants/tmp_in_gene.bed ]; then
        awk 'BEGIN{OFS="\t"} {
            chrom  = $1
            pos    = $3    # END = original 1-based POS
            gene   = $7
            print gene, chrom, pos
        }' variants/tmp_in_gene.bed \
        | sort -k1,1 -k2,2 -k3,3n \
        | uniq \
        > variants/tmp_in_gene_key.txt

        # Join with original VCF info using CHROM+POS as key
        # Create lookup from vcf_info: key=CHROM_POS value=REF ALT QUAL FILTER
        awk 'BEGIN{OFS="\t"} {
            key = $1"_"$2
            val = $3"\t"$4"\t"$5"\t"$6
            lookup[key] = val
        }
        END {
            # This block intentionally empty
        }' variants/tmp_vcf_info.txt > /dev/null

        # Better approach: process both files together with awk
        awk 'BEGIN{OFS="\t"}
        NR==FNR {
            # Reading vcf_info: CHROM POS REF ALT QUAL FILTER
            key = $1"_"$2
            info[key] = $3"\t"$4"\t"$5"\t"$6
            next
        }
        {
            # Reading in_gene.bed: CHROM START END CHROM_B START_B END_B GeneID . Strand
            chrom = $1
            pos   = $3
            gene  = $7
            key   = chrom"_"pos
            if (key in info) {
                print gene, chrom, pos, info[key]
            }
        }' variants/tmp_vcf_info.txt variants/tmp_in_gene.bed \
        | sort -u \
        >> "$TSV_OUT"
    fi

    # Process intergenic variants
    if [ -s variants/tmp_intergenic.bed ]; then
        awk 'BEGIN{OFS="\t"}
        NR==FNR {
            # Reading vcf_info
            key = $1"_"$2
            info[key] = $3"\t"$4"\t"$5"\t"$6
            next
        }
        {
            # Reading intergenic.bed: CHROM START END
            chrom = $1
            pos   = $3
            key   = chrom"_"pos
            if (key in info) {
                print "intergenic", chrom, pos, info[key]
            }
        }' variants/tmp_vcf_info.txt variants/tmp_intergenic.bed \
        | sort -u \
        >> "$TSV_OUT"
    fi

    # Cleanup
    rm -f variants/tmp_vcf_info.txt \
          variants/tmp_vcf_sorted.bed \
          variants/tmp_gene_sorted.bed \
          variants/tmp_in_gene.bed \
          variants/tmp_intergenic.bed \
          variants/tmp_in_gene_key.txt \
          2>/dev/null

    # Final stats
    local TOTAL IN_GENE INTERGENIC
    TOTAL=$(tail -n +2 "$TSV_OUT" | wc -l)
    IN_GENE=$(tail -n +2 "$TSV_OUT" | awk '$1 != "intergenic"' | wc -l)
    INTERGENIC=$(tail -n +2 "$TSV_OUT" | awk '$1 == "intergenic"' | wc -l)

    echo "  --- Final Stats ---"
    echo "  Total annotated rows:     $TOTAL"
    echo "  Variants in gene regions: $IN_GENE"
    echo "  Intergenic variants:      $INTERGENIC"
    if [ "$TOTAL" -gt 0 ]; then
        awk -v ig="$IN_GENE" -v t="$TOTAL" \
            'BEGIN{printf "  Gene region coverage:     %.1f%%\n", ig/t*100}'
    fi
    echo "  -------------------"
    echo ""
    echo "  Output first 8 lines:"
    head -8 "$TSV_OUT"
    echo ""
}

# --- Step 3: Annotate raw VCF ---
echo "[2/4] Annotating raw VCF (all.raw.vcf.gz)..."
annotate_vcf \
    "variants/all.raw.vcf.gz" \
    "variants/all.raw.annotated_geneid.tsv" \
    "Raw VCF"

# --- Step 4: Annotate filtered VCF ---
echo "[3/4] Annotating filtered VCF (all.filtered.vcf.gz)..."
annotate_vcf \
    "variants/all.filtered.vcf.gz" \
    "variants/all.filtered.annotated_geneid.tsv" \
    "Filtered VCF"

echo "[4/4] All done!"
echo ""
echo "====================================================="
echo "  Output files:"
echo "  variants/all.raw.annotated_geneid.tsv"
echo "  variants/all.filtered.annotated_geneid.tsv"
echo "====================================================="