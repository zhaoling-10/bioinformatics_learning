#!/usr/bin/env bash
set -euo pipefail

module load samtools 2>/dev/null || true
module load picard 2>/dev/null || true

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
LE_DICT="${LE_REF%.fa}.dict"
LT_DICT="${LT_REF%.fa}.dict"

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

ensure_ref_support_files() {
    local ref_fa="$1"
    local ref_dict="$2"

    if [[ ! -f "${ref_fa}.fai" ]]; then
        echo "Index missing for ${ref_fa}; creating .fai with samtools faidx"
        samtools faidx "${ref_fa}"
    fi

    if [[ ! -f "${ref_dict}" ]]; then
        echo "Dictionary missing for ${ref_fa}; creating ${ref_dict} with Picard"
        picard CreateSequenceDictionary \
            R="${ref_fa}" \
            O="${ref_dict}"
    fi
}

ensure_ref_support_files "${LE_REF}" "${LE_DICT}"
ensure_ref_support_files "${LT_REF}" "${LT_DICT}"

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

if [[ ! -f "${CHAINDIR}/LE_lifted_to_LT.vcf.gz" || ! -f "${CHAINDIR}/LT_lifted_to_LE.vcf.gz" ]]; then
    echo "ERROR: Liftover did not produce expected output VCF files."
    echo "Check logs:"
    echo "  - ${CHAINDIR}/liftover_LE_to_LT.log"
    echo "  - ${CHAINDIR}/liftover_LT_to_LE.log"
    exit 1
fi

# Check results
echo "=== LE→LT liftover results ==="
if command -v bcftools >/dev/null 2>&1; then
    bcftools stats "${CHAINDIR}/LE_lifted_to_LT.vcf.gz" | grep "^SN" || true
    echo "=== Number of sites that failed to liftover ==="
    bcftools view -H "${CHAINDIR}/LE_to_LT_rejected.vcf.gz" | wc -l
else
    echo "WARNING: bcftools not found, skip stats summary."
    echo "Install/load bcftools and run manually:"
    echo "  bcftools stats ${CHAINDIR}/LE_lifted_to_LT.vcf.gz | grep '^SN'"
    echo "  bcftools view -H ${CHAINDIR}/LE_to_LT_rejected.vcf.gz | wc -l"
fi