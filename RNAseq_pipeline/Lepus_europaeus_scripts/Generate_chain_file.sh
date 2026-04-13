# ============================================================
# Set variables (modify according to your actual file names)
# ============================================================
WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LE_REF="${WORKDIR}/lepus_europaeus_gff/ref/genome.fa"   # LE reference genome
LT_REF="${WORKDIR}/lepus_timidus_gff/ref/genome.fa"     # LT reference genome

CHAINDIR="${WORKDIR}/liftover"

mkdir -p ${CHAINDIR}
cd ${CHAINDIR}

module load minimap2 2>/dev/null || true
module load samtools 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true

# ============================================================
# Method: Use minimap2 to do whole genome alignment, then convert to chain file
# Need to install: minimap2, samtools, ucsc-tools(or use Python script to convert)
# ============================================================

# Step 2a：Use minimap2 to do LE_REF → LT_REF whole genome alignment
# -cx asm5: suitable for same-species or closely related species alignment (sequence similarity >95%)
# -cx asm10: suitable for species with 90-95% similarity (more appropriate for these two Lepus species)

minimap2 \
    -cx asm10 \
    --cs \
    -t 10 \
    ${LT_REF} \
    ${LE_REF} \
    > LE_to_LT.paf

echo "minimap2 alignment completed"

# Step 2b: Convert PAF format to a chain file
# Use paftools (bundled with minimap2)
paftools.js chain LE_to_LT.paf > LE_to_LT.chain

# If paftools.js is unavailable, use this Python script instead:
python3 paf_to_chain.py LE_to_LT.paf LE_to_LT.chain

echo "Chain file generation completed"

# Similarly, generate the reverse chain (LT → LE)
minimap2 \
    -cx asm10 \
    --cs \
    -t 10 \
    ${LE_REF} \
    ${LT_REF} \
    > LT_to_LE.paf

paftools.js chain LT_to_LE.paf > LT_to_LE.chain