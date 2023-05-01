#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=00:30:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Slurm_files/Step2.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Slurm_files/Step2.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_scripts/2_spearman_corrs.R
