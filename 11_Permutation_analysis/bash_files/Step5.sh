#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=2:00:00
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Slurm_files/Step5.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Slurm_files/Step5.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_scripts/5_Plot_dists_genes_clusters_overlaps.R
