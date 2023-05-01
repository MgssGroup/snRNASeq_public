#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step3.5.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Slurm_files/Step3.5.out

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_scripts/3.5_readd_stability_info.R
