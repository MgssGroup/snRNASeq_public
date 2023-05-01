#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH --time=2:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step2.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step2.out

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_scripts/2_cluster_annotation.R
