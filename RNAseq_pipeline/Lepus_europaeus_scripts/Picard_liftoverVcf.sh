#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Base paths on CSC server
# ============================================================
WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
CHAINDIR="${WORKDIR}/liftover"
LE_REF="${WORKDIR}/lepus_europaeus_gff/ref/genome.fa"
LT_REF="${WORKDIR}/lepus_timidus_gff/ref/genome.fa"

# ============================================================
# Confirm VCF file paths
# ============================================================
LE_VCF="${WORKDIR}/lepus_europaeus_gff/variants/all.filtered.vcf.gz"
LT_VCF="${WORKDIR}/lepus_timidus_gff/variants/all.filtered.vcf.gz"

module load picard 2>/dev/null || true
module load bcftools 2>/dev/null || true

if [[ ! -f "${LE_VCF}" || ! -f "${LT_VCF}" ]]; then
    echo "ERROR: input VCF not found."
    echo "Expected:"
    echo "  - ${LE_VCF}"
    echo "  - ${LT_VCF}"
    exit 1
fi

if [[ ! -f "${CHAINDIR}/LE_to_LT.chain" || ! -f "${CHAINDIR}/LT_to_LE.chain" ]]; then
    echo "ERROR: chain file not found."
    echo "Expected:"
    echo "  - ${CHAINDIR}/LE_to_LT.chain"
    echo "  - ${CHAINDIR}/LT_to_LE.chain"
    exit 1
fi

# If your VCF names are different, inspect available files first.
ls "${WORKDIR}/lepus_europaeus_gff/variants/" "${WORKDIR}/lepus_timidus_gff/variants/"

# ============================================================
# Method A: Picard LiftoverVcf (recommended)
# ============================================================

# LE_ref → LT_ref
picard LiftoverVcf \
    INPUT="${LE_VCF}" \
    OUTPUT="${CHAINDIR}/LE_lifted_to_LT.vcf.gz" \
    CHAIN="${CHAINDIR}/LE_to_LT.chain" \
    REJECT="${CHAINDIR}/LE_to_LT_rejected.vcf.gz" \
    REFERENCE_SEQUENCE="${LT_REF}" \
    RECOVER_SWAPPED_REF_ALT=true \
    WRITE_ORIGINAL_POSITION=true \
    2>&1 | tee "${CHAINDIR}/liftover_LE_to_LT.log"

# LT_ref → LE_ref
picard LiftoverVcf \
    INPUT="${LT_VCF}" \
    OUTPUT="${CHAINDIR}/LT_lifted_to_LE.vcf.gz" \
    CHAIN="${CHAINDIR}/LT_to_LE.chain" \
    REJECT="${CHAINDIR}/LT_to_LE_rejected.vcf.gz" \
    REFERENCE_SEQUENCE="${LE_REF}" \
    RECOVER_SWAPPED_REF_ALT=true \
    WRITE_ORIGINAL_POSITION=true \
    2>&1 | tee "${CHAINDIR}/liftover_LT_to_LE.log"

# Check results
echo "=== LE→LT liftover results ==="
bcftools stats "${CHAINDIR}/LE_lifted_to_LT.vcf.gz" | grep "^SN"
echo "=== Number of sites that failed to liftover ==="
bcftools view -H "${CHAINDIR}/LE_to_LT_rejected.vcf.gz" | wc -l