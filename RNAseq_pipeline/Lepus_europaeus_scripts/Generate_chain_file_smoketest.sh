#!/usr/bin/env bash
set -euo pipefail

# Smoke test for chain generation:
# 1) Build small FASTA subsets from the original genomes
# 2) Run minimap2 alignment in both directions
# 3) Convert PAF to chain
# 4) Print quick QC summary

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LE_REF="${WORKDIR}/lepus_europaeus_gff/ref/genome.fa"
LT_REF="${WORKDIR}/lepus_timidus_gff/ref/genome.fa"

TESTDIR="${WORKDIR}/liftover_smoketest"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Keep conservative defaults for close species.
MM2_THREADS="${MM2_THREADS:-4}"
MM2_INDEX_SIZE="${MM2_INDEX_SIZE:-2G}"
MM2_PRESET="${MM2_PRESET:-asm10}"

# Subset mode:
# - contigs: keep first N_CONTIGS records
# - bp: keep records until total sequence length reaches TARGET_BP
SUBSET_MODE="${SUBSET_MODE:-bp}"
N_CONTIGS="${N_CONTIGS:-10}"
TARGET_BP="${TARGET_BP:-50000000}"

mkdir -p "${TESTDIR}/logs"
cd "${TESTDIR}"
RUN_LOG="logs/smoketest_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

module load minimap2 2>/dev/null || true
module load samtools 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true

echo "Smoke test started at: $(date)"
echo "Runtime log: ${TESTDIR}/${RUN_LOG}"
echo "SUBSET_MODE=${SUBSET_MODE}, N_CONTIGS=${N_CONTIGS}, TARGET_BP=${TARGET_BP}"
echo "MM2_PRESET=${MM2_PRESET}, MM2_THREADS=${MM2_THREADS}, MM2_INDEX_SIZE=${MM2_INDEX_SIZE}"

subset_fasta() {
    local in_fa="$1"
    local out_fa="$2"
    local n="$3"
    awk -v n="${n}" '
        /^>/ {h++; if (h > n) exit}
        h > 0 && h <= n {print}
    ' "${in_fa}" > "${out_fa}"
}

subset_fasta_by_bp() {
    local in_fa="$1"
    local out_fa="$2"
    local target_bp="$3"
    awk -v target="${target_bp}" '
        BEGIN {bp=0; keep=0}
        /^>/ {
            if (bp >= target) exit
            keep=1
            print
            next
        }
        {
            if (keep) {
                print
                bp += length($0)
            }
        }
    ' "${in_fa}" > "${out_fa}"
}

convert_paf_to_chain() {
    local paf_file="$1"
    local chain_file="$2"
    if command -v paftools.js >/dev/null 2>&1 && paftools.js 2>&1 | grep -q -w "chain"; then
        paftools.js chain "${paf_file}" > "${chain_file}"
    else
        python3 "${SCRIPT_DIR}/paf_to_chain.py" "${paf_file}" "${chain_file}"
    fi
}

if [[ "${SUBSET_MODE}" == "bp" ]]; then
    subset_fasta_by_bp "${LE_REF}" "LE_subset.fa" "${TARGET_BP}"
    subset_fasta_by_bp "${LT_REF}" "LT_subset.fa" "${TARGET_BP}"
    echo "Subset FASTA files created by total base pairs"
elif [[ "${SUBSET_MODE}" == "contigs" ]]; then
    subset_fasta "${LE_REF}" "LE_subset.fa" "${N_CONTIGS}"
    subset_fasta "${LT_REF}" "LT_subset.fa" "${N_CONTIGS}"
    echo "Subset FASTA files created by contig count"
else
    echo "ERROR: SUBSET_MODE must be 'bp' or 'contigs'"
    exit 1
fi

minimap2 -x "${MM2_PRESET}" -c -t "${MM2_THREADS}" -I "${MM2_INDEX_SIZE}" \
    "LT_subset.fa" "LE_subset.fa" > LE_to_LT.test.paf
convert_paf_to_chain LE_to_LT.test.paf LE_to_LT.test.chain

minimap2 -x "${MM2_PRESET}" -c -t "${MM2_THREADS}" -I "${MM2_INDEX_SIZE}" \
    "LE_subset.fa" "LT_subset.fa" > LT_to_LE.test.paf
convert_paf_to_chain LT_to_LE.test.paf LT_to_LE.test.chain

echo "Smoke test outputs:"
ls -lh LE_subset.fa LT_subset.fa LE_to_LT.test.paf LE_to_LT.test.chain LT_to_LE.test.paf LT_to_LE.test.chain
echo "PAF line counts:"
wc -l LE_to_LT.test.paf LT_to_LE.test.paf
echo "First 2 PAF lines (LE -> LT):"
awk 'NR<=2{print}' LE_to_LT.test.paf

echo "Smoke test finished at: $(date)"
echo "Runtime log saved to: ${TESTDIR}/${RUN_LOG}"
