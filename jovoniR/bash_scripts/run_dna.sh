#!/bin/bash
#SBATCH --job-name=devil_DNA_pdos
#SBATCH --partition=GENOA
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --output=logs/devil_DNA.out
#SBATCH --error=logs/devil_DNA.err
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process

# 2. Run the R script
Rscript fit_DNA.R