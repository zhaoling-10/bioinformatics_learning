Instructions for these scripts:

1. 00_setup.sh --> Project initialization and directory creation
2. 01_config_sratookit.sh --> Configuring SRA tools
3. 02_prefetch.sh --> Downloading raw data from NCBI
4. 03_fasterq_dump.sh --> SRA to FASTQ format conversion
5. 04_prepare_reference.sh --> Preparing reference genome and annotations
6. 05_fastq.sh --> Quality control, adapter removal, and pruning
7. 06_star_index.sh --> Building the STAR alignment index
8. 07_star_align.sh --> Sequence alignment to the genome
9. 08_post_alignment_qc.sh --> Alignment quality statistics
10. 09_featurecounts.sh --> Gene expression counting
11. 10_variant_calling.sh --> Variance detection (this script takes a long time)
12. 11_multiqc.sh --> Generating a comprehensive quality report

Because I couldn't use the prefetch function on CSC, so I downloaded data in fastq format from EBI. The script name downloading the data is download_with_enaDataGet_smart.sbatch, and I didn't try running the 01-03 script.
