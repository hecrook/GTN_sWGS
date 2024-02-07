#!/bin/bash
#SBATCH --job-name=index-nextflow
#SBATCH --partition=master-worker
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=hannah.crook@icr.ac.uk
#SBATCH --mail-type=ALL

srun nextflow 6-GTN_sWGS.nf -resume