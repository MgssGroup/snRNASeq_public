#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --gres=gpu:k80:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --error=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step4.err
#SBATCH --output=/home/malosree/scratch/Finalized_merging_normalization/Slurm_files/Step4.out

Rscript /home/malosree/scratch/Finalized_merging_normalization/Finalized_scripts/4_basic_plots.R
