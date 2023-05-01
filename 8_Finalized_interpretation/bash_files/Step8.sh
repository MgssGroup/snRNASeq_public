#!/bin/bash
#SBATCH --account=def-gturecki
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step8.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step8.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_scripts/8_per_subject_per_cluster_stats.R
