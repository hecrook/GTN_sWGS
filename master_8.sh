#!/bin/bash
#SBATCH --job-name=index-nextflow
#SBATCH --partition=master-worker
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=hannah.crook@icr.ac.uk
#SBATCH --mail-type=ALL

srun nextflow 8-GTN_sWGS.nf -resume --binsize 500 --genome hg38
