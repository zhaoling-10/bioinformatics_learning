#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bcftools 2>/dev/null || true

cd "$WORKDIR"

# one step:Generate site information using bcftools mpileup and pass it to bcftools call via a pipeline for variant detection.
bcftools mpileup \
  -Ou \
  -f ref/genome.fa \
  -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  aln/*.Aligned.sortedByCoord.out.bam \
| bcftools call -mv -Ou \
| bcftools view -Oz -o variants/all.raw.vcf.gz
# Ou --> Output uncompressed BCF format (faster)
# f ref/genome.fa --> Reference genome FASTA file
# q 20 -Q 20 --> Set minimum mapping quality (-q) and minimum base quality (-Q) to 20
# a DP,DP,AD --> Add total depth (DP) and allele depth (AD) information
# aln/*.Aligned.sortedByCoord.out.bam --> All alignments generate BAM files
# bcftools call -mv -Ou --> Call variants (-m: multi-allelic mode, -v: output only variant sites)
# bcftools view -Oz -o variants/all.raw.vcf.gz --> Convert to compressed VCF format

# create orginal VCF files index
bcftools index -f variants/all.raw.vcf.gz

# two step: Filter variants (Quality >= 30, Depth >= 5)
bcftools filter \
  -i 'QUAL>=30 && INFO/DP>=5' \
  -Oz -o variants/all.filtered.vcf.gz \
  variants/all.raw.vcf.gz
# Filter conditions
# Output the filtered compressed VCF
# Enter the original VCF file

# Create an index for the filtered VCF file
bcftools index -f variants/all.filtered.vcf.gz

# Output the number of filtered variants
echo -n "Filtered variant count: "
bcftools view -H variants/all.filtered.vcf.gz | wc -l