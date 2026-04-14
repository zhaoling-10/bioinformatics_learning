#!/usr/bin/env bash
#SBATCH --job-name=generate_chain
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

PROJECT_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
# Script location on the server (same project directory as input data).
GENERATE_CHAIN_SCRIPT="${PROJECT_DIR}/Generate_chain_file.sh"

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

module load minimap2 2>/dev/null || true
module load samtools 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true

# Publication/final-analysis profile (conservative defaults for close species):
# - Keep asm10 preset
# - Keep extra minimap2 options empty (do not filter secondary hits by default)
export MM2_THREADS=8
export MM2_INDEX_SIZE=8G
export MM2_PRESET=asm10
export MM2_EXTRA_OPTS=""

echo "Job started at: $(date)"
echo "Running publication profile:"
echo "  MM2_PRESET=${MM2_PRESET}"
echo "  MM2_EXTRA_OPTS='${MM2_EXTRA_OPTS}'"
echo "  MM2_THREADS=${MM2_THREADS}"
echo "  MM2_INDEX_SIZE=${MM2_INDEX_SIZE}"
echo "Slurm stdout log: ${PROJECT_DIR}/logs/generate_chain_${SLURM_JOB_ID}.out"
echo "Slurm stderr log: ${PROJECT_DIR}/logs/generate_chain_${SLURM_JOB_ID}.err"

bash "${GENERATE_CHAIN_SCRIPT}"

echo "Job finished at: $(date)"
echo "Check logs:"
echo "  - ${PROJECT_DIR}/logs/generate_chain_${SLURM_JOB_ID}.out"
echo "  - ${PROJECT_DIR}/logs/generate_chain_${SLURM_JOB_ID}.err"
