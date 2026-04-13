# ============================================================
# Set variables (modify according to your actual file names)
# ============================================================
#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LE_REF="${WORKDIR}/lepus_europaeus_gff/ref/genome.fa"   # LE reference genome
LT_REF="${WORKDIR}/lepus_timidus_gff/ref/genome.fa"     # LT reference genome

CHAINDIR="${WORKDIR}/liftover"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Allow overriding to reduce memory usage when needed
MM2_THREADS="${MM2_THREADS:-4}"
MM2_INDEX_SIZE="${MM2_INDEX_SIZE:-4G}"

mkdir -p "${CHAINDIR}"
cd "${CHAINDIR}"

module load minimap2 2>/dev/null || true
module load samtools 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true

# ============================================================
# Method: Use minimap2 to do whole genome alignment, then convert to chain file
# Need to install: minimap2, samtools, ucsc-tools(or use Python script to convert)
# ============================================================

convert_paf_to_chain() {
    local paf_file="$1"
    local chain_file="$2"

    if command -v paftools.js >/dev/null 2>&1 && paftools.js 2>&1 | grep -q -w "chain"; then
        paftools.js chain "${paf_file}" > "${chain_file}"
        echo "Converted ${paf_file} to ${chain_file} using paftools.js"
    else
        echo "paftools.js chain is unavailable, using Python fallback converter"
        python3 "${SCRIPT_DIR}/paf_to_chain.py" "${paf_file}" "${chain_file}"
    fi
}

# Step 2a：Use minimap2 to do LE_REF → LT_REF whole genome alignment
# -cx asm5: suitable for same-species or closely related species alignment (sequence similarity >95%)
# -cx asm10: suitable for species with 90-95% similarity (more appropriate for these two Lepus species)

minimap2 \
    -cx asm10 \
    --cs \
    -t "${MM2_THREADS}" \
    -I "${MM2_INDEX_SIZE}" \
    "${LT_REF}" \
    "${LE_REF}" \
    > LE_to_LT.paf

echo "minimap2 alignment completed"

# Step 2b: Convert PAF format to a chain file
# Use paftools (bundled with minimap2)
convert_paf_to_chain LE_to_LT.paf LE_to_LT.chain

echo "Chain file generation completed"

# Similarly, generate the reverse chain (LT → LE)
minimap2 \
    -cx asm10 \
    --cs \
    -t "${MM2_THREADS}" \
    -I "${MM2_INDEX_SIZE}" \
    "${LE_REF}" \
    "${LT_REF}" \
    > LT_to_LE.paf

convert_paf_to_chain LT_to_LE.paf LT_to_LE.chain