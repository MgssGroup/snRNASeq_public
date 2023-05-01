#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step1.3.1.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step1.3.1.out

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_scripts/1.3.1_morabito_proportions.R
