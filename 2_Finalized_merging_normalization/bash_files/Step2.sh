#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --gres=gpu:k80:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=30G
#SBATCH --time=03:00:00
#SBATCH --error=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step2.err
#SBATCH --output=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step2.out

Rscript /home/malosree/scratch/Finalized_merging_normalization/Finalized_scripts/2_merge_normalize.R
