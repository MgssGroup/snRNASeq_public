#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=1:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step4.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step4.out

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_scripts/4_cluster_markers_plots.R
