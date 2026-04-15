#!/usr/bin/env bash
#SBATCH --job-name=generate_chain_prod
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

PROJECT_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
SCRIPT_DIR="/home/zhaoling/bioinformatics_learning/RNAseq_pipeline/Lepus_europaeus_scripts"
CHAIN_SCRIPT="${SCRIPT_DIR}/Generate_chain_file_production.sh"
CHECK_SCRIPT="${SCRIPT_DIR}/check_liftover_toolchain_csc.sh"

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

if [[ ! -f "${CHAIN_SCRIPT}" ]]; then
    echo "ERROR: production chain script not found: ${CHAIN_SCRIPT}"
    exit 1
fi

if [[ ! -f "${CHECK_SCRIPT}" ]]; then
    echo "ERROR: toolchain check script not found: ${CHECK_SCRIPT}"
    exit 1
fi

module load minimap2 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true
module load samtools 2>/dev/null || true

# Publication/production profile defaults
export MM2_PRESET="${MM2_PRESET:-asm10}"
export MM2_THREADS="${MM2_THREADS:-${SLURM_CPUS_PER_TASK}}"
export MM2_INDEX_SIZE="${MM2_INDEX_SIZE:-8G}"

echo "Job started at: $(date)"
echo "Running preflight toolchain check..."
bash "${CHECK_SCRIPT}"
echo "Running production chain generation..."
echo "MM2_PRESET=${MM2_PRESET} MM2_THREADS=${MM2_THREADS} MM2_INDEX_SIZE=${MM2_INDEX_SIZE}"
echo "Slurm stdout log: ${PROJECT_DIR}/logs/generate_chain_prod_${SLURM_JOB_ID}.out"
echo "Slurm stderr log: ${PROJECT_DIR}/logs/generate_chain_prod_${SLURM_JOB_ID}.err"

bash "${CHAIN_SCRIPT}"

echo "Job finished at: $(date)"
