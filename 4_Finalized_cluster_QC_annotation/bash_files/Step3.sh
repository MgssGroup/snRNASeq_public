#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step3.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step3.out

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_scripts/3_cluster_matching.R
