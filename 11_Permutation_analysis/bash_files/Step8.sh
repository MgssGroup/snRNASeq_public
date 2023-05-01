#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --time=3:00:00
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Slurm_files/Step8.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Slurm_files/Step8.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_scripts/8_spearman_corrs_permuted.R
