#!/bin/bash
#SBATCH --job-name=index-nextflow
#SBATCH --partition=master-worker
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=hannah.crook@icr.ac.uk
#SBATCH --mail-type=ALL

srun nextflow main.nf --binsize 500 --genome hg38 --input GTN2022_sWGS_samplesheet_batch02.csv
