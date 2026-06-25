#!/bin/bash
#SBATCH --job-name=index-nextflow
#SBATCH --partition=master-worker
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=hannah.crook@icr.ac.uk
#SBATCH --mail-type=ALL

srun nextflow main.nf \
--binsize 500 \
--genome hg38 \
--input GTN2022_sWGS_samplesheet_batch02.csv \
--fasta /data/reference-data/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
--adap /data/scratch/DMP/UCEC/EVOLIMMU/hcrook/gestational_sWGS/nextflow/bin/TruSeq3-PE-2.fa \
--bedfile /data/scratch/DMP/UCEC/EVOLIMMU/hcrook/gestational_sWGS/nextflow/docs/ETT31-sWGS11_sorted_genome_windows_1kb.bed \
-resume
