#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=2:00:00
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Slurm_files/Step16.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Slurm_files/Step16.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_scripts/16_pseudotime_replots_source_data.R
