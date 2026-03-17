#!/usr/bin/env bash
set -euo pipefail   # Security Settings:-e: Immediately exit if any command in the script fails;-u: Report an error when using an undefined variable;-o pipefail: If any command in a pipeline fails, the entire pipeline is considered to have failed

# ====== EDIT ONLY THESE IF NEEDED ======
WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/lepus_europaeus_gff"
THREADS=10

# Your reference files (must be present in $WORKDIR/ref_src/)
REF_FA_GZ="GCF_033115175.1_mLepTim1.pri_genomic.fna.gz"
REF_GTF_GZ="Lepus_timidus_annotation.zip"
# REF_FA_GZ="GCA_040893245.2_mLepTim1.1_pri_genomic.fna.gz"
# REF_GTF_GZ="GCA_040893245.2_mLepTim1.1_pri_genomic.gtf.gz"

# STAR sjdbOverhang (readLength-1; your data are 2x101bp -> 100)
SJDB_OVERHANG=100
# =======================================

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Project structure
mkdir -p ref_src ref sra fastq trimmed star_index aln qc counts variants logs tmp

# Create list.txt if not present
if [ ! -f list.txt ]; then
cat > list.txt <<'EOF'
SRR18740835
SRR18740836
SRR18740837
SRR18740838
SRR18740839
SRR18740840
SRR18740841
SRR18740842
EOF
fi

# Save config file used by all steps
cat > config.sh <<EOF
WORKDIR="${WORKDIR}"
THREADS=${THREADS}
REF_FA_GZ="${REF_FA_GZ}"
REF_GTF_GZ="${REF_GTF_GZ}"
SJDB_OVERHANG=${SJDB_OVERHANG}
EOF

echo "Setup done in: $WORKDIR"
echo "Next:"
echo "1) Copy your reference files into: $WORKDIR/ref_src/"
echo "   - $REF_FA_GZ"
echo "   - $REF_GTF_GZ"
echo "2) Run: bash 01_config_sratoolkit.sh"