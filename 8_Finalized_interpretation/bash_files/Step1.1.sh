#!/bin/bash
#SBATCH --account=def-gturecki
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step1.1.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step1.1.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_scripts/1.1_fgsea.R
