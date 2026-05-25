#!/bin/bash
#SBATCH --job-name=RNA_boot_pdos
#SBATCH --partition=GENOA
#SBATCH --mem=200GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --output=logs/devil_RNA_boot.out
#SBATCH --error=logs/devil_RNA_boot.err
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process
Rscript /orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/diploid_bootstrap_RNA.R