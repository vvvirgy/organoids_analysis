#!/bin/bash
#SBATCH --job-name=devil_organoid
#SBATCH --partition=GENOA
#SBATCH --mem=200GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/devil_RNA.out
#SBATCH --error=logs/devil_RNA.err

source ~/.bashrc
conda activate process

# Run the R script
Rscript /orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/fit_RNA.R psinorm FALSE