#!/usr/bin/env bash
set -euo pipefail

# Production chain generation script.
# NOTE: This script intentionally does NOT use paf_to_chain.py fallback.
# It only runs when paftools.js has a real 'chain' subcommand.

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LE_REF="${WORKDIR}/lepus_europaeus_gff/ref/genome.fa"
LT_REF="${WORKDIR}/lepus_timidus_gff/ref/genome.fa"
CHAINDIR="${WORKDIR}/liftover_prod"

MM2_THREADS="${MM2_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
MM2_INDEX_SIZE="${MM2_INDEX_SIZE:-8G}"
MM2_PRESET="${MM2_PRESET:-asm10}"

mkdir -p "${CHAINDIR}/logs"
cd "${CHAINDIR}"
RUN_LOG="logs/generate_chain_prod_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

module load minimap2 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true
module load samtools 2>/dev/null || true

echo "Run started at: $(date)"
echo "Runtime log: ${CHAINDIR}/${RUN_LOG}"
echo "MM2_PRESET=${MM2_PRESET} MM2_THREADS=${MM2_THREADS} MM2_INDEX_SIZE=${MM2_INDEX_SIZE}"

if ! command -v paftools.js >/dev/null 2>&1; then
    echo "ERROR: paftools.js not found."
    echo "Run check script first:"
    echo "  bash check_liftover_toolchain_csc.sh"
    exit 1
fi

if ! paftools.js 2>&1 | grep -q -w "chain"; then
    echo "ERROR: paftools.js exists but does not support 'chain' subcommand."
    echo "This production script refuses Python fallback to avoid incompatible chain output."
    echo "Please switch to a toolchain where paftools.js chain is available."
    exit 1
fi

echo "Confirmed mature chain path: minimap2 + paftools.js chain"

minimap2 -x "${MM2_PRESET}" -c -t "${MM2_THREADS}" -I "${MM2_INDEX_SIZE}" \
    "${LT_REF}" "${LE_REF}" > LE_to_LT.paf
paftools.js chain LE_to_LT.paf > LE_to_LT.chain
echo "LE -> LT done"

minimap2 -x "${MM2_PRESET}" -c -t "${MM2_THREADS}" -I "${MM2_INDEX_SIZE}" \
    "${LE_REF}" "${LT_REF}" > LT_to_LE.paf
paftools.js chain LT_to_LE.paf > LT_to_LE.chain
echo "LT -> LE done"

echo "Output summary:"
ls -lh LE_to_LT.paf LE_to_LT.chain LT_to_LE.paf LT_to_LE.chain
echo "Run finished at: $(date)"
echo "Runtime log saved to: ${CHAINDIR}/${RUN_LOG}"

