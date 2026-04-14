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
# Allow overriding for performance/memory trade-offs
MM2_THREADS="${MM2_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
MM2_INDEX_SIZE="${MM2_INDEX_SIZE:-4G}"
MM2_PRESET="${MM2_PRESET:-asm10}"
# Keep defaults conservative for closely related species.
MM2_EXTRA_OPTS="${MM2_EXTRA_OPTS:-}"
read -r -a MM2_EXTRA_OPTS_ARR <<< "${MM2_EXTRA_OPTS}"

mkdir -p "${CHAINDIR}"
cd "${CHAINDIR}"
mkdir -p logs
RUN_LOG="logs/generate_chain_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1
echo "Run started at: $(date)"
echo "Runtime log: ${CHAINDIR}/${RUN_LOG}"
echo "minimap2 preset/options: -x ${MM2_PRESET} ${MM2_EXTRA_OPTS}"
echo "minimap2 threads/index: -t ${MM2_THREADS} -I ${MM2_INDEX_SIZE}"

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
    -x "${MM2_PRESET}" \
    -c \
    "${MM2_EXTRA_OPTS_ARR[@]}" \
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
    -x "${MM2_PRESET}" \
    -c \
    "${MM2_EXTRA_OPTS_ARR[@]}" \
    -t "${MM2_THREADS}" \
    -I "${MM2_INDEX_SIZE}" \
    "${LE_REF}" \
    "${LT_REF}" \
    > LT_to_LE.paf

convert_paf_to_chain LT_to_LE.paf LT_to_LE.chain
echo "Run finished at: $(date)"
echo "Runtime log saved to: ${CHAINDIR}/${RUN_LOG}"