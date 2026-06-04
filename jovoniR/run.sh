#!/usr/bin/env bash
#SBATCH --job-name=
#SBATCH --partition=GENOA
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

source ~/.bashrc
conda activate process

Rscript 01_fit_data.R
