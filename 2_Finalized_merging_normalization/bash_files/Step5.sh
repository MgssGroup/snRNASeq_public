#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --gres=gpu:k80:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#SBATCH --time=4:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Slurm_files/Step5.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Slurm_files/Step5.out

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_scripts/5_reharmonize_seeded.R
