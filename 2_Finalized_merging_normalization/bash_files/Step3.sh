#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --gres=gpu:k80:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --time=04:00:00
#SBATCH --error=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step3.err
#SBATCH --output=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step3.out

Rscript /home/malosree/scratch/Finalized_merging_normalization/Finalized_scripts/3_harmonize.R
