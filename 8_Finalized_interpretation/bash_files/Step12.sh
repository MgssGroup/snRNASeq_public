#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step12.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Slurm_files/Step12.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_scripts/12_CellChat_Mic_PV.R
