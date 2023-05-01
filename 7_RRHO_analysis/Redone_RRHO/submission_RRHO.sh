#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Redone_RRHO/submission_RRHO.out
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Redone_RRHO/submission_RRHO.err

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Redone_RRHO/run_RROH2_per_clusterID_sample_2022.04.17.R
