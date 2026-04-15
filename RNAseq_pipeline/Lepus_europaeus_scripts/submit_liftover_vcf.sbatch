#!/usr/bin/env bash
#SBATCH --job-name=liftover_vcf
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

PROJECT_DIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LIFTOVER_SCRIPT="${PROJECT_DIR}/Picard_liftoverVcf.sh"

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

if [[ ! -f "${LIFTOVER_SCRIPT}" ]]; then
    echo "ERROR: liftover script not found: ${LIFTOVER_SCRIPT}"
    exit 1
fi

module load picard 2>/dev/null || true
module load samtools 2>/dev/null || true
if ! command -v bcftools >/dev/null 2>&1; then
    module load samtools 2>/dev/null || true
fi

# JVM settings for compute node execution.
export _JAVA_OPTIONS="${_JAVA_OPTIONS:--Xms1g -Xmx8g}"

echo "Job started at: $(date)"
echo "Running script: ${LIFTOVER_SCRIPT}"
echo "Slurm stdout log: ${PROJECT_DIR}/logs/liftover_vcf_${SLURM_JOB_ID}.out"
echo "Slurm stderr log: ${PROJECT_DIR}/logs/liftover_vcf_${SLURM_JOB_ID}.err"

bash "${LIFTOVER_SCRIPT}"

echo "Job finished at: $(date)"
echo "Check logs:"
echo "  - ${PROJECT_DIR}/logs/liftover_vcf_${SLURM_JOB_ID}.out"
echo "  - ${PROJECT_DIR}/logs/liftover_vcf_${SLURM_JOB_ID}.err"