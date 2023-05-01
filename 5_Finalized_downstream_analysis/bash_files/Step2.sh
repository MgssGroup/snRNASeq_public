#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=5:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step2.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step2.out

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_scripts/2_pseudotime_comparison.R
