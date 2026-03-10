#!/bin/bash
#SBATCH --job-name=devil_array
#SBATCH --partition=GENOA
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --output=logs/mutlimodal_%A_%a.out
#SBATCH --error=logs/mutlimodal_%A_%a.err
#SBATCH --array=1-4
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process

# 1. Get the parameters for this specific task
PARAMS_FILE="params/RNA_config.txt"
METHOD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -d' ' -f1)
STABLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -d' ' -f2)

echo "Task ID $SLURM_ARRAY_TASK_ID: Running $METHOD with use_stable=$STABLE"

# 2. Run the R script
Rscript analysis_cheng_et_al.R "$METHOD" "$STABLE"